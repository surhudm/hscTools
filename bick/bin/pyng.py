#!/usr/bin/env python
# Original filename: pyng.py
#
# Author: 
# Email: 
# Date: Thu 2014-07-31 18:38:03
# 
# Summary: 
# 

import sys
import os
import re
import math
import argparse

import multiprocessing

import subprocess
import time
import datetime
import signal


def colorize(text, color, bold=False):

    colors = {
        "red"    :"31",
        "green"  :"32",
        "yellow" :"33",
        "blue"   :"34",
        "magenta":"35",
        "cyan"   :"36",
        "grey"   :"37",
        }
    
    base = "\033["
    code = colors[color]
    if bold:
        code += ";1"

    prefix = base + code + "m"
    suffix = base + "0m"
    return prefix + text + suffix


class Color(object):
    def __init__(self, active=False):
        self.active = active
        
    def __call__(self, fmt, value):
        if self.active:
            if value > 0.0:
                c = "green"
            else:
                c = "red"
            return colorize(fmt%value, c)
        else:
            return fmt%value

        

def all_done():
    print "Received SIGINT.  Exiting."
    sys.exit()

def init():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def ping(node):
    try:
        output = subprocess.check_output(['ping', '-c', '1', '-W', '1', node], stderr=subprocess.STDOUT)
    except:
        output = ""
        
    t = 0.0
    m = re.search("time=(\d+\.\d+)", output)
    if m:
        t = float(m.groups()[0])
    return t


#############################################################
#
# Main body of code
#
#############################################################

def main(network, n, startn=1, alive=False, color=False, line=False, interval=1, timestamp=False):

    nodes = ["%s.%d"%(network, x) for x in range(startn,startn+n)]

    c_obj = Color(color)
    
    pool_size = multiprocessing.cpu_count() * 2
    if len(nodes) < pool_size:
        pool_size = len(nodes)

    fmt = "%03d" if alive else "%5d"
    h = [fmt%(x) for x in range(startn,startn+n)]

    marg = 0
    if timestamp:
        marg = 9
    
    if alive:
        heads = []
        for i in range(3):
            heads.append(" "*marg + " ".join([x[i] for x in h]))
        header = "\n".join(heads)
        s_len = len(heads[0])
    else:
        header = " "*marg + " ".join(h)
        s_len = len(header)

    print header
    print '-'*s_len
        
        
    while(True):

        pool = multiprocessing.Pool(pool_size, init)
        try:
            pool_outputs = pool.map(ping, nodes)
        except KeyboardInterrupt:
            pool.terminate()
            pool.join()
            all_done()
        else:
            pool.close() # no more tasks
            pool.join()

        t = ""
        if timestamp:
            t = datetime.datetime.now().strftime("%H:%M:%S ")
        try:
            if alive:
                s = [c_obj("%d",math.ceil(x)) for x in pool_outputs]
            else:
                s = [c_obj("%5.3f",x) for x in pool_outputs]
            print t+" ".join(s)
            if line:
                print "-------------"
            time.sleep(interval)
        except KeyboardInterrupt:
            all_done()
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action='version', version="0.0")
    parser.add_argument("network", help="IP net e.g. 192.168.100")
    parser.add_argument("n", type=int, help="Number of hosts")
    parser.add_argument("-i", "--interval", type=int, default=1, help="Interval between pings in seconds.")
    parser.add_argument("-s", "--startn", type=int, default=1, help="Number to start at")
    parser.add_argument("-a", "--alive", default=False, action='store_true',
                        help="Just show 1/0 (deal/alive)")
    parser.add_argument("-c", '--color', default=False, action='store_true', help="Use color")
    parser.add_argument("-l", "--line", default=False, action='store_true', help="Separate lines with ---")
    parser.add_argument("-t", "--timestamp", default=False, action='store_true', help="Print a timestamp")
    args = parser.parse_args()

    main(args.network, args.n, args.startn,
         alive=args.alive, color=args.color, line=args.line, interval=args.interval, timestamp=args.timestamp)
