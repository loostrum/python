import argparse

argparser=argparse.ArgumentParser(description='Process some ints.')
argparser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
argparser.add_argument('--sum', dest='accumulate',action='store_const', const=sum, default=max, help='sum the integers (default:find the max)')
args=argparser.parse_args()
print args.accumulate()

