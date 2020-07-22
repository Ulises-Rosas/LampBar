import argparse
import re

parser = argparse.ArgumentParser(description= "", add_help=False)

parser.add_argument('-l', metavar='lamp',
                    help='LAMP primers from laval software')

parser.add_argument('-p', metavar='csv',
                    help='position in a csv file')

args = parser.parse_args()

primers = open(args.l, "r").readlines()
pos     = open(args.p, "r").readline().replace("\n","")


whole_matches = []

for i in primers[::16]:

    tmp  = i.replace("\n","")
    tmp2 = re.sub( ".*\(locations: (.*)\)","\\1", tmp )
    tmp3 = re.findall("[0-9]+-[0-9]+", tmp2)
    # tmp3 = ['262-284', '468-489', '302-321', '446-466', '366-388', '392-412']

    p_matches = []

    for tmp4 in tmp3:
        ## [262, 284]
        tmp5 = [ int(y) for y in  tmp4.split('-')]

        matches = 0

        for tmp_pos in pos.split(","):

            if int(tmp_pos) in range( tmp5[0], tmp5[1] ):

                matches += 1

        p_matches.append( matches )

    whole_matches.append( p_matches )

sum_wp = []
sum_p  = []

for tmp_sum in [ list( reversed(tmp_rev) ) for tmp_rev in whole_matches ]:

    sum_wp.append( str( sum( tmp_sum ) ) )

    sum_p.append(
        ":#match per primer: " + "(" + ",".join( [ str(l) for l in tmp_sum ] ) + ")"
    )

iterThat = zip(primers[9::16 ],
               primers[11::16],
               primers[5::16 ],
               primers[7::16 ],
               primers[1::16 ],
               primers[3::16 ],
               sum_wp,
               sum_p )

for f1,b1,f2,b2,f3,b3,ws,s in iterThat:
    print(
        ",".join(

            [f1.replace("\n",""),
             b1.replace("\n",""),
             f2.replace("\n",""),
             b2.replace("\n",""),
             f3.replace("\n",""),
             b3.replace("\n",""),
             ws,
             s]
        )
    )
