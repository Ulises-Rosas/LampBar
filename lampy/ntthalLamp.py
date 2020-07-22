import subprocess
import os
import argparse
import itertools
import re


path = subprocess.Popen(['which', 'ntthal'],
                        stdout = subprocess.PIPE,
                        stderr = subprocess.STDOUT).stdout.read().decode('utf-8').replace("\n","")

# help_path
if path is "":
    help_path = 'path to ntthal program'
else:
    help_path = 'path to ntthal program. Default = [%s]' % path


parser = argparse.ArgumentParser(description = """

            -----------------------------------------------------------------------------
            ============== Simple wrapper of ntthal program for LAMP primers ============
            -----------------------------------------------------------------------------
            
                                
            This script combinatorially evaluates the presence of both hairpins and dimers under 
            below parameters. 
                                """, add_help = True)

parser.add_argument('-l',
                    metavar = 'primers',
                    help    = """
                    List of lamp primers in csv format.
                    Structure of input: F1,R1,F2,R2,F3,R3.
                    """)
parser.add_argument('-p',
                    metavar = 'path',
                    help    = help_path,
                    default = path)
parser.add_argument('-f',
                    action  = 'store_true',
                    help    = 'if selected, this transforms six-primer format into three-primer format (i.e. B3,F3 BIP, FIP)'
                    )
parser.add_argument('-e',
                    action  = 'store_true',
                    help    = 'if selected, print out explanation of filtering'
                    )

parser.add_argument('-ht',
                    metavar = 'cut-off',
                    help    = 'Filter Option: Primer hairpin Tm threshold...[default = 40]',
                    default = 40)
parser.add_argument('-Hd',
                    metavar = 'cut-off',
                    help    = 'Filter Option: Primer hairpin dG threshold...[default = -10000]',
                    default = -10000)
parser.add_argument('-dt',
                    metavar = 'cut-off',
                    help    = 'Filter Option: Cross dimer Tm threshold......[default = 40]',
                    default = 40)
parser.add_argument('-Dd',
                    metavar = 'cut-off',
                    help    = 'Filter Option: Cross dimer dG threshold......[default = -10000]',
                    default = -10000)


parser.add_argument('-mv',
                    metavar = 'mM',
                    help    = 'ntthal Option: Monovalent cations in mM......[default = 50]',
                    default = 50)
parser.add_argument('-dv',
                    metavar = 'mM',
                    help    = 'ntthal Option: Divalent cations in mM........[default = 0]',
                    default = 0)
parser.add_argument('-n',
                    metavar = 'mM',
                    help    = 'ntthal Option: concentration of dNTP in mM...[default = 0]',
                    default = 0)
parser.add_argument('-c',
                    metavar = 'nM',
                    help    = 'ntthal Option: concentration of dNTP in mM...[default = 50]',
                    default = 50)

parser.add_argument('-t',
                    metavar = 'Tm ',
                    help    = 'ntthal Option: Temperature for calculations..[default = 62.5 C]',
                    default = 62.5)

args = parser.parse_args()

if args.p is None or args.p is "" :
    print("\nNo path of ntthal found. Try to find it with: \033[1;35mwhich ntthal\033[0m")
    os._exit(1)

if args.l is None:
    print("\nPlease provide a list of primer")
    os._exit(1)

def revcom(string):
    """
    BIP: revcom(R1) + R2
    FIP: revcom(F1) + F2
    """
    libco = {"A": "T", "G": "C", "C": "G",
             "T": "A", "R": "Y", "Y": "R",
             "S": "S", "W": "W", "K": "M",
             "M": "K", "B": "V", "V": "B",
             "D": "H", "H": "D", "N": "N",
             "-": "-"}

    complement = []
    for pos in range(0, len(string)):
        complement.append(libco[string[pos]])

    return "".join(complement)[::-1]


path = str(args.p)
file = open(str(args.l), "r").readlines()

## arbitrary
# path = "/Users/admin/primer3-2.4.0/src/ntthal"
# file = open("listOfLamps", "r").readlines()
## arbitrary

ntthal_args = [path,
               "-mv", str(args.mv),
               "-dv", str(args.dv),
               "-n" , str(args.n),
               "-d" , str(args.c),
               "-t" , str(args.t),
               "-s1"]

convert   = args.f
explainMe = args.e

for lamp in file:

    ## delete
    # lamp = file[0]
    ## delete

    # print("\n\ttesting: %s\n" % lamp)

    p_tmp = lamp.replace("\n", "").split(",")

    ## heres where information is loaded
    # p_tmp[6:]
    ## heres where information is loaded

    libLamp = {'f3'      : p_tmp[4],
               'b3'      : p_tmp[5],
               'bip_tmp' : revcom( p_tmp[1] ) + p_tmp[3],
               'fip_tmp' : revcom( p_tmp[0] ) + p_tmp[2]}



    explanations = []

    ### hairpis
    for h in itertools.combinations( libLamp.values(), 1):


        complete_args = ntthal_args + [ h[0] ] + ["-a", "HAIRPIN"]

        # print(" ".join(complete_args))

        result = subprocess.Popen(complete_args,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE).stdout.read().decode('utf-8')

        dG = re.findall("dG = [-0-9.]+",result)
        t  = re.findall("t = [-0-9.]+" ,result)

        if dG.__len__() > 0 and t.__len__() > 0:

            dG = float( re.sub("dG = ([-0-9.]+)", "\\1", dG[0]) )
            t =  float( re.sub("t = ([-0-9.]+)" , "\\1",  t[0]) )

            if dG < float(args.Hd) and t <= float(args.ht):
                h_issue = [a for a, b in libLamp.items() if re.findall(b, h[0])][0]
                explanations.append("dG of %s < %s (%s)" % (h_issue,float(args.Hd), dG))

            elif dG >= float(args.Hd) and t > float(args.ht):
                h_issue = [a for a, b in libLamp.items() if re.findall(b, h[0])][0]
                explanations.append("Tm of %s > %s (%s)" % (h_issue, float(args.ht), t))

            elif dG < float(args.Hd) and t > float(args.ht):

                h_issue = [a for a, b in libLamp.items() if re.findall(b, h[0])][0]
                explanations.append("dG of %s < %s (%s) and Tm of %s > %s (%s)" % (h_issue, float(args.Hd), dG,
                                                                                       h_issue, float(args.ht), t))

    for d in itertools.combinations(libLamp.values(), 2):

        complete_args =  ntthal_args + [d[0], "-s2", d[1]]

        # print(" ".join(complete_args))

        result = subprocess.Popen(complete_args,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE).stdout.read().decode('utf-8')

        complete_args1 = ntthal_args + [d[0], "-s2", d[1], "-a", "END1"]

        # print(" ".join(complete_args1))

        result1 = subprocess.Popen(complete_args1,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE).stdout.read().decode('utf-8')

        complete_args2 = ntthal_args + [d[0], "-s2", d[1], "-a", "END2"]

        # print(" ".join(complete_args2))

        result2 = subprocess.Popen(complete_args2,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE).stdout.read().decode('utf-8')

        super_results = result + result1 + result2

        dG = re.findall("dG = [-0-9.]+", super_results)
        t  = re.findall("t = [-0-9.]+" , super_results)

        if dG.__len__() > 0 and t.__len__() > 0:

            allTm = [ float( re.sub("t = ([-0-9.]+)",  "\\1", i) )  for i in t ]
            # print(allTm)
            alldG = [ float( re.sub("dG = ([-0-9.]+)", "\\1", i) )  for i in dG]
            # print(alldG)


            anyGreater_Tm = any( i > float(args.dt) for i in allTm )
            anyLesser_dG  = any( i < float(args.Dd) for i in alldG )


            if anyLesser_dG and not anyGreater_Tm:

                d_issue = [a for a,b in libLamp.items() if re.findall("(%s)" % "|".join(d), b)]
                dG      = ",".join([i for i in alldG if i < float(args.Dd)])

                explanations.append("dG of %s + %s < %s (%s)" % (d_issue[0], d_issue[1], float(args.Dd), dG))

            elif not anyLesser_dG and anyGreater_Tm:

                d_issue = [a for a, b in libLamp.items() if re.findall("(%s)" % "|".join(d), b)]
                Tm = ",".join([i for i in allTm if i > float(args.dt) ])

                explanations.append("Tm of %s + %s > %s (%s)" % (d_issue[0], d_issue[1], float(args.dt), Tm))

            elif anyLesser_dG and anyGreater_Tm:

                d_issue = [a for a, b in libLamp.items() if re.findall("(%s)" % "|".join(d), b)]

                Tm = ",".join([i for i in allTm if i > float(args.dt) ])
                dG = ",".join([i for i in alldG if i < float(args.Dd) ])

                explanations.append("dG of %s + %s < %s (%s) and Tm of %s + %s > %s (%s)" % (d_issue[0], d_issue[1], float(args.Dd), dG,
                                                                                             d_issue[0], d_issue[1], float(args.dt), Tm))

    tmp_ex = ",".join(explanations)

    if explainMe and convert:

        if tmp_ex == "":

            print(",".join(libLamp.values()) + ", :" +  ",".join(libLamp.keys()) + ": ," + ",".join(p_tmp[6:]) + ", Ok" )
        else:
            print(",".join(libLamp.values()) + ", :" +  ",".join(libLamp.keys()) + ": ," + ",".join(p_tmp[6:]) +  ", " + tmp_ex )

    elif not explainMe and convert:
        if tmp_ex == "":
            print(",".join(libLamp.values()) + ", :" +  ",".join(libLamp.keys()) + ": ," + ",".join(p_tmp[6:]))

    elif explainMe and not convert:
        if tmp_ex == "":
            print(lamp.replace("\n", "") + ", Ok" )
        else:
            print(lamp.replace("\n", "") + ", " + tmp_ex)

    elif not explainMe and not convert:
        if tmp_ex == "":
            print( lamp.replace("\n", "") )

