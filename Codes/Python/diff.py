#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/10/04 16:16
# @Author   : Iydon
# @File     : course4.1.py

import re

class diff:
    def __init__(self, variable:list):
        self.variable = variable
        self.diff_map = {"sin(_)": "cos(_)",
                         "cos(_)": "-sin(_)",
                         "tan(_)": "1/cos(_)**2",
                         "exp(_)": "exp(_)",
                         "_**n"  : "n*_**(n-1)",
                         "log(_)": "1/_"}

    def diff(self, formula, target:str=""):
        self.target = target if target else self.variable[0]

    def test(self):
        print(self.__split_once("1*2","*"))


    def __diff_element(self, expr):
        return self.diff_map(expr)

    def __split_once(self, expr, sep):
        seper = "(?<=[^\{0}])\{0}(?=[^\{0}])"
        return re.split(seper.format(sep), expr)

    def __split_time(self, expr):
        return re.split("(?<=[^\*])\*(?=[^\*])", expr)

    def __split_plus(self, expr):
        return re.split("(?<=[^\+])\+(?=[^\+])", expr)

    def __target_variable(self, expr):
        prepost = "[^.\d]"
        pattern = "(?<={0}){1}(?={0})".format(prepost, self.target)
        return re.sub(pattern, "_", expr)

d = diff(["x"])
d.diff("sin(x)")
d.test()