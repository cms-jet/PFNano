import os
import argparse

"""
Parse and return the arguments provided by the user.
"""

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--command", dest="command", required=True,
    choices=["status", 'resubmit', 'kill'],
    help=("The command to run over crab dirs"))
parser.add_argument("-d", "--dir", dest="crab_dir", required=True,
    help=("The crab directory you want to use"))

parser.add_argument("--blacklist", dest="blacklist", type=str,
    help=("Extension to append"))

parser.add_argument("--maxjobruntime", dest="maxjobruntime", type=str,
    help=("Extra time"))

parser.add_argument("--run", dest="run", action='store_true',
    help=("Execute"))

args = parser.parse_args()

for d in sorted(os.listdir(args.crab_dir)):
    line = "crab {} {}/{}".format(args.command, args.crab_dir, d)
    if args.blacklist is not None:
        line += " --siteblacklist " + args.blacklist
    if args.maxjobruntime is not None:
        line += " --maxjobruntime " + args.maxjobruntime
    if args.command == 'status':
        print(d)
        line += " | grep -e Jobs -e finished -e idle -e running -e transferring -e 'CRAB server:'| grep -v Warning"
    if args.run:
        os.system(line)
    else:
        print(line)

if not args.run:
    print("\nPass '--run' to execute the above commands.")

    
