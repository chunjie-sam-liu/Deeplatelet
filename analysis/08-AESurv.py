# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-01-13 20:30:01
# @DESCRIPTION:

import numpy as np
import feather
import matplotlib.pyplot as plt

# For preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn_pandas import DataFrameMapper

# torch modules
import torch
from torch import nn
import torch.nn.functional as F
from torch.nn.modules.linear import Linear
import torchtuples as tt

# pycox modules
from pycox.models import LogisticHazard, PMF, DeepHit, MTLR, BCESurv
from pycox.models.loss import NLLLogistiHazardLoss
from pycox.evaluation import EvalSurv

import os

root_dir = os.path.dirname(os.path.abspath(__file__))

# random
np.random.seed(1234)
_ = torch.manual_seed(1234)


class NetAESurv(nn.Module):
    def __init__(self, in_features, encoded_features, out_features):
        super().__init__()

        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(in_features, 2048),
            nn.ReLU(),
            nn.Linear(2048, 1024),
            nn.ReLU(),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, encoded_features),
        )

        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(encoded_features, 128),
            nn.ReLU(),
            nn.Linear(128, 256),
            nn.ReLU(),
            nn.Linear(256, 512),
            nn.ReLU(),
            nn.Linear(512, 1024),
            nn.ReLU(),
            nn.Linear(1024, 2048),
            nn.ReLU(),
            nn.Linear(2048, in_features),
        )

        # Full connection
        self.survnet = nn.Sequential(
            nn.Linear(encoded_features, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Linear(32, out_features),
        )

    def forward(self, input):
        encoded = self.encoder(input)
        decoded = self.decoder(encoded)
        phi = self.survnet(encoded)
        return phi, decoded

    def predict(self, input):
        encoded = self.encoder(input)
        return self.survnet(encoded)


class LossAELogHaz(nn.Module):
    def __init__(self, alpha):
        super().__init__()
        assert (alpha >= 0) and (alpha <= 1), "Need `alpha` in [0, 1]."
        self.alpha = alpha
        self.loss_surv = NLLLogistiHazardLoss()
        self.loss_ae = nn.MSELoss()

    def forward(self, phi, decoded, target_loghaz, target_ae):
        idx_durations, events = target_loghaz
        loss_surv = self.loss_surv(phi, idx_durations, events)
        loss_ae = self.loss_ae(decoded, target_ae)
        return self.alpha * loss_surv + (1 - self.alpha) * loss_ae


def load_data(filepath):
    df = feather.read_dataframe(source=filepath)
    df_train = df.loc[df.oc == "OC521"].drop(columns=["barcode", "oc"], axis=1)
    df_val = df.loc[df.oc == "OC44"].drop(columns=["barcode", "oc"], axis=1)
    df_test1 = df.loc[df.oc == "OC79"].drop(columns=["barcode", "oc"], axis=1)
    df_test2 = df.loc[df.oc == "OC172"].drop(columns=["barcode", "oc"], axis=1)
    return df_train, df_val, df_test1, df_test2


def get_target(df):
    return (df["duration"].values, df["event"].values)


def transform_features(df_train, df_val, df_test1, df_test2):
    columns = df_train.columns
    columns = columns[: len(columns) - 2]
    standardize = [([col], StandardScaler()) for col in columns]

    x_mapper = DataFrameMapper(standardize)

    x_train = x_mapper.fit_transform(df_train).astype("float32")
    x_val = x_mapper.transform(df_val).astype("float32")
    x_test1 = x_mapper.transform(df_test1).astype("float32")
    x_test2 = x_mapper.transform(df_test2).astype("float32")

    return x_train, x_val, x_test1, x_test2


def transform_labels(df_train, df_val, nd=10):
    num_durations = nd
    labtrans = LogisticHazard.label_transform(num_durations)
    y_train_surv = labtrans.fit_transform(*get_target(df_train))
    y_val_surv = labtrans.transform(*get_target(df_val))

    return y_train_surv, y_val_surv, labtrans


def run(filepath):
    # load data
    df_train, df_val, df_test1, df_test2 = load_data(filepath)
    # transform features
    x_train, x_val, x_test1, x_test2 = transform_features(df_train, df_val, df_test1, df_test2)
    # transform labels
    y_train_surv, y_val_surv, labtrans = transform_labels(df_train, df_val)

    # make train and validation datasets with tuplefy
    train = tt.tuplefy(x_train, (y_train_surv, x_train))
    val = tt.tuplefy(x_val, (y_val_surv, x_val))

    duration_test1, events_test1 = get_target(df_test1)
    duration_test2, events_test2 = get_target(df_test2)

    # set arch
    in_features = x_train.shape[1]
    encoded_features = 64
    out_features = labtrans.out_features
    netaesurv = NetAESurv(in_features, encoded_features, out_features)
    print(netaesurv)

    # loss
    loss = LossAELogHaz(0.6)
    print(loss)

    # model
    model = LogisticHazard(net=netaesurv, optimizer=tt.optim.Adam(0.01), duration_index=labtrans.cuts, loss=loss)

    # dl = model.make_dataloader(train, batch_size=5, shuffle=False)
    # batch = next(iter(dl))
    # model.compute_metrics(batch)

    # metrics
    metrics = dict(loss_surv=LossAELogHaz(1), loss_ae=LossAELogHaz(0))

    # callbacks
    callbacks = [tt.cb.EarlyStopping()]

    # cycling
    batch_size = 5
    epochs = 3

    # trainning model
    log = model.fit(
        *train, batch_size=batch_size, epochs=epochs, callbacks=callbacks, verbose=True, val_data=val, metrics=metrics
    )

    res = model.log.to_pandas()
    _ = res[["train_loss", "val_loss"]].plot()


def main():
    run(filepath="data/rda/total416.os.se.norm.coxph.feather")


if __name__ == "__main__":
    main()
