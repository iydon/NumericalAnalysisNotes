#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time     : 2018/09/19 09:00
# @Author   : Iydon
# @File     : Iydonize.py
import re


def str2pipe(string:str, identifier:str="|>", first_flag:bool=True):
    """
    管道操作符.
    ----------------------------
    Args:
        identifier: Default, `|>'.
        first_flag: Bool, 1|f(x)=>f(1,x).

    Returns:
        Determinted by the last function.

    Raises:
        None.
    """
    def __pipe(*target):
        """
        Tuple -> object.
        """
        result = target[0][0](*target[0][-1])
        if len(target)==1:
            return result
        else:
            if first_flag:
                last = target[1][0],(result,)+target[1][-1]
            else:
                last = target[1][0],target[1][-1]+(result,)
            return __pipe(*((last,)+target[2:]))

    def __extract_function(string:str) -> str:
        """
        提取函数名.
        """
        if "(" in string and ")" in string:
            return string[:string.index("(")]
        else:
            return string

    def __extract_args(string:str) -> str:
        """
        提取参数.
        """
        if "(" in string and ")" in string:
            idx = lambda s: string.index(s)
            return "(%s,)"%string[idx("(")+1:idx(")")]
        else:
            return "()"

    __str = re.sub(" +"," ",string).split(identifier)
    e_f,e_a = __extract_function,__extract_args
    result = ["(%s,%s)"%(e_f(__s),e_a(__s)) for __s in __str]
    return eval("__pipe(%s)"%(",".join(result)))



result = """
max(1,2,3) |> lambda x: x**2 |> print
"""
str2pipe(string=result)
