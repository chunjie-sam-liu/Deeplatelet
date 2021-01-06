# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-01-13 20:30:01
# @DESCRIPTION:

import numpy as np
import feather

# For preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn_pandas import DataFrameMapper

import torch
from torch import nn
import torch.nn.functional as F
import torchtuples as tt

from pycox.models import LogisticHazard
from pycox.models.loss import NLLLogistiHazardLoss
from pycox.evaluation import EvalSurv

import os

root_dir = os.path.dirname(os.path.abspath(__file__))

# random
np.random.seed(1234)
_ = torch.manual_seed(1234)


def run():
    pass


def main():
    run()


if __name__ == "__main__":
    main()
