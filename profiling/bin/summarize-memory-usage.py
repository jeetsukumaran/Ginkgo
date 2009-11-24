#! /usr/bin/env python

import os
import sys
import re
from datetime import datetime
from optparse import OptionGroup
from optparse import OptionParser

def extract_memory_usage(fpath, rss=True, subsample=1):
    results = {}
    for idx, row in enumerate(open(fpath, "rU")):
        if idx == 0 or subsample==0 or ((idx % subsample) != 0):
            continue
        items = [i for i in re.split('\s+', row) if i]
        if len(items) != 8:
            # incomplete row
            break
        elapsed_time = items[3].split(':')
        t = float(elapsed_time[0]) + float(elapsed_time[1])/60
        if len(elapsed_time) > 2:
            t += float(elapsed_time[2])/3600
        rss = items[6]
        vsize = items[7]
        if rss:
            results[t] = rss
        else:
            results[t] = vsize
    return results

def main():
    """
    Main CLI handler.
    """

    parser = OptionParser(add_help_option=True, usage="%prog [options] <GINKGO OUTPUT LOGS>")

    parser.add_option('-r', '--rss',
        action='store_true',
        dest='rss',
        default=True,
        help='extract Resident State Set (i.e., physical) memory usage (default)')

    parser.add_option('-v', '--vsize',
        action='store_true',
        dest='rss',
        default=False,
        help='extract virtual memory memory usage')

    parser.add_option('-s', '--subsample',
        action='store',
        dest='subsample',
        default=600,
        help='subsample step (default=%default)')

    parser.add_option('-u', '--unstack',
        action='store_false',
        dest='stack',
        default=True,
        help='unstack output (separate column for every population size)')

    parser.add_option('-b', '--break-on-na',
        action='store_true',
        dest='break_on_na',
        default=False,
        help='terminate output when first NA is encountered (default=%default)')

    (opts, args) = parser.parse_args()

    psizes = []
    results = {}
    for a in args:
        sys.stderr.write("Processing '%s' ...\n" % a)
        psize = os.path.basename(a).split('.')[0].split('_')[0]
        psizes.append(psize)
        results[psize] = extract_memory_usage(a, opts.rss, opts.subsample)

    memory_sizes = []
    for r in results.values():
        memory_sizes.extend(r.keys())
    memory_sizes = list(set(memory_sizes))
    memory_sizes.sort()
    sys.stdout.write('Gen\t')
    if opts.stack:
        sys.stdout.write('Memory\tPopSize\n')
    else:
        sys.stdout.write('%s\n' % ('\t'.join(psizes)))
    na_found = False
    for memsz in memory_sizes:
        if opts.break_on_na and na_found:
            break
        row = []
        if opts.stack:
            for p in psizes:
                r = results[p]
                if memsz in r:
                    sys.stdout.write('%s\t%s\t%s\n' % (memsz, r[memsz], p[1:]))
                else:
                    if opts.break_on_na:
                        na_found = True
                        break
        else:
            for p in psizes:
                r = results[p]
                if memsz in r:
                    row.append(str(r[memsz]))
                else:
                    if opts.break_on_na:
                        na_found = True
                        break
                    else:
                        row.append('NA')
            if not (opts.break_on_na and na_found):
                sys.stdout.write('%s\t' % memsz)
                sys.stdout.write('%s\n' % ('\t'.join(row)))

if __name__ == '__main__':
    main()


