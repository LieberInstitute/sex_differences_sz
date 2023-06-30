## Run each command in parallel

import argparse
import pandas as pd
from os import system,environ

def get_cmd(task, tissue):
    return pd.read_csv(f"{tissue}.assemble.cmd",
                       header=None).iloc[task-1,0]


def run_sge(task, tissue):
    cmd = get_cmd(task, tissue)
    system(cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=int, required=True, metavar="N",
                        help="The task number. Assumes starts with 1.")
    parser.add_argument('--tissue', type=str, default="caudate",
                        help="Tissue to run [default: %(default)s].")
    args=parser.parse_args()
    run_sge(args.task, args.tissue)


if __name__ == '__main__':
    main()
