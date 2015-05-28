#!/usr/bin/env python

import sys, os, re
import argparse
import subprocess
import collections

Job = collections.namedtuple("Job", "user queue jobname sessid nds tsk memreq timreq status elaptime")

class Record(object):
    def __init__(self, name, jobs=None):
        self.name = name
        self.data = {"np": 0, "jobs": {}}
        self.jobs = jobs
        
    def update(self, line):
        line = line.strip()
        if len(line) == 0:
            return

        fields = line.split(" = ")
        if len(fields) == 2:
            k, v = fields
            if k == 'jobs':
                clean = [re.sub("^\s*\d+\/", "", x) for x in v.split(',')]
                v = collections.Counter(clean)
            if k in 'np':
                v = int(v)
            self.data[k] = v

    def __str__(self):
        node_str = ""
        jobs = 0
        for jobid,count in self.data['jobs'].items():
            if self.jobs:
                node_str += self.jobs[jobid].user.lower()[0]*count
            else:
                node_str += "*"*count
            jobs += count
        
        node = node_str + "-" * (self.data['np'] - jobs)
        s = "%8s %-16s %02d/%02d  %s" % (self.name, self.data['state'], self.data['np'], jobs, node)
        return s

def main():
    queue = subprocess.check_output(("qstat", "-a"), stderr=subprocess.STDOUT)
    jobs = dict()
    for line in queue.split("\n"):
        if re.match("^\d", line):
            fields = line.split()
            jobs[fields[0]] = Job(*fields[1:])

            
    output = subprocess.check_output("pbsnodes", stderr=subprocess.STDOUT)

    record = None
    records = []
    for line in output.split('\n'):
        if re.match("^\w", line):
            record = Record(line.strip(), jobs)
            records.append(record)
        else:
            record.update(line)

    for r in sorted(records, key=lambda x: x.name):
        print r

    print ""
    count = collections.defaultdict(int)
    for jid, j in jobs.items():
        count[j.user] += 1

    abbreviations = set()
    for k,v in count.items():
        abbrev = k[0]
        if abbrev in abbreviations:
            abbrev = k[0].upper()
        print "%1s %-16s: %d" % (abbrev, "("+k+")", v)
        abbreviations.add(abbrev)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    main()
    
