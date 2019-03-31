import argparse
import re
import os
import operator

parser = argparse.ArgumentParser(description= "", add_help=False)

parser.add_argument('-a', metavar='Seqs',
                    help='Sequences')

parser.add_argument('-p', metavar='csv',
                    help='position in a csv file')

parser.add_argument('-r',
                    action = 'store_true',
                    help='''Reverse analysis. That is, from no hyphen position to allocate
                    diagnostic nucleotides into a new alignment''')

args = parser.parse_args()

def fas_to_dic(x):
    ##fasta file is read
    file = open(x)

    ##all lines are stored in a variable
    file_content = file.readlines()

    ##an empty variable is created to store lines without the
    ##newline metacharacter symbol (i.e. "\n")
    seqs_list = []

    ##Given this following loop:
    for i in file_content:
        ##newline symbol is
        ##replaced by nothing
        seqs_list.append(i.replace("\n", ""))

        ##here is where keys are going to be stored
    keys = []

    ##likewise, here is where values are going to be stored
    values = []

    i = 0
    ##first value for dictionary's keys
    while (">" in seqs_list[i]):

        keys.append(seqs_list[i])
        i += 1
        JustOneValue = []

        while ((">" in seqs_list[i]) == False):
            JustOneValue.append(seqs_list[i])
            i += 1

            if (i == len(seqs_list)):
                i -= 1
                break
        values.append("".join(JustOneValue))

    return dict(zip(keys, values))

complete_aln = fas_to_dic(str(args.a))

pos     = "|".join( [ "^%s$" % i for i in open( str(args.p) ).readline().replace("\n","").split(",") ] )

lengths = { i:list(j.replace("-","")).__len__() for i,j in complete_aln.items() }
s_item  = sorted( lengths.items(), key = operator.itemgetter(1) )[-1][0]


if args.r:

    old_pos = 0
    new_pos = []

    for p,q in enumerate(complete_aln[s_item]):
        if q != "-":
            old_pos += 1
            if re.findall( "(" + pos + ")" , str(old_pos)):
                new_pos.append(p + 1)

else:
    new_pos = []
    tmp_pos = 0

    for p, q in enumerate(complete_aln[s_item]):
        if q != "-":
            tmp_pos += 1
            if re.findall("(" + pos + ")", str(p + 1)):
                new_pos.append(tmp_pos)


print(
    ",".join([str(i) for i in new_pos])
)
