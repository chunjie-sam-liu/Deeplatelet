# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-01-18 11:49:37
# @DESCRIPTION:


def debug(func):
    def wrapper(*args, **kwargs):
        print("Prepare and say")
        return func(*args, **kwargs)

    return wrapper


@debug
def say(something):
    print(something)


say("hi")
