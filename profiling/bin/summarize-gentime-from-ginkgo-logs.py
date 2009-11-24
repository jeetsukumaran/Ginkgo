#! /usr/bin/env python

import os
import sys
import re
from datetime import datetime
from optparse import OptionGroup
from optparse import OptionParser

def extract_generation_times(fpath):
    pattern = re.compile('\[(.*)\] Generation (\d+) life-cycle running.')
    rows = open(fpath, "rU").readlines()
    results = {}
    start_time = None
    for row in rows:
        m = pattern.match(row)
        if (m):
            timestamp_str = m.groups(1)[0]
            timestamp = datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M:%S")
            generation = int(m.groups(1)[1])
            if generation == 0:
                start_time = timestamp
                results[0] = 0
            else:
                hours, mins, secs = str(timestamp - start_time).split(":")
                results[generation] = float(hours) + float(mins)/60 + float(secs)/3600
    return results

def main():
    """
    Main CLI handler.
    """

    parser = OptionParser(add_help_option=True, usage="%prog [options] <GINKGO OUTPUT LOGS>")

    parser.add_option('--stop-on-na',
        action='store_true',
        dest='stop_on_na',
        default=False,
        help='terminate output when first NA is encountered (default=%default)')

    (opts, args) = parser.parse_args()

    psizes = []
    results = {}
    for a in args:
        sys.stderr.write("Processing '%s' ...\n" % a)
        psize = os.path.basename(a).split('.')[0].split('_')[0]
        psizes.append(psize)
        results[psize] = extract_generation_times(a)

    gens = []
    for r in results.values():
        gens.extend(r.keys())
    gens = list(set(gens))
    gens.sort()
    sys.stdout.write('Gen\t')
    sys.stdout.write('%s\n' % ('\t'.join(psizes)))
    na_found = False
    for g in gens:
        if opts.stop_on_na and na_found:
            break
        row = []
        for p in psizes:
            r = results[p]
            if g in r:
                row.append(str(r[g]))
            else:
                if opts.stop_on_na:
                    na_found = True
                    break
                else:
                    row.append('NA')
        if not (opts.stop_on_na and na_found):
            sys.stdout.write('%s\t' % g)
            sys.stdout.write('%s\n' % ('\t'.join(row)))

if __name__ == '__main__':
    main()


