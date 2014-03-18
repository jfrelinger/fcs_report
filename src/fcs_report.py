'''
Created on Mar 18, 2014

@author: Jacob Frelinger <jacob.frelinger@gmail.com>
'''

import fcm
from glob import glob
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt


def set_color(bp):
    ''' set boxplot colors'''
    plt.setp(bp['boxes'][0], color='b')
    plt.setp(bp['caps'][0], color='b')
    plt.setp(bp['caps'][1], color='b')
    plt.setp(bp['whiskers'][0], color='b')
    plt.setp(bp['whiskers'][1], color='b')
    plt.setp(bp['fliers'][0], color='b')
    plt.setp(bp['fliers'][1], color='b')
    plt.setp(bp['medians'][0], color='r')

def hsh(x):
    '''create a semi-unique id out of a fcmdata object'''
    return hash(''.join(x.long_names))


def report(glb, out):
    '''generate a nice summary report of a group of fcs files'''
    panels = {}
    panel_map = {}
    means = defaultdict(list)
    nevents = defaultdict(list)
    laserabs = defaultdict(list)
    for i in glob(glb):
        print i
        x = fcm.loadFCS(i)
        h = hsh(x)
        if h not in panels:
            panels[h] = x.long_names
        for i in range(len(x.channels)):
            laserabs[x.short_names[i]].append(x.channels[i])

        panel_map[i] = h
        means[h].append(x.mean(0))
        nevents[h].append(x.shape[0])

    panel_idx = {}
    for i, j in enumerate(panels.keys()):
        panel_idx[j] = i

    rev_map = defaultdict(list)
    for i in panel_map:
        rev_map[panel_idx[panel_map[i]]].append(i)

    common = set(panels[panels.keys()[0]])
    for i in panels:
        common.intersection_update(set(panels[i]))


    with open('panels.md', 'w') as f:
        f.write('Panels\n')
        f.write('======\n')
        for i, j in enumerate(panels):
            f.write('Panel %d:\n' % i)
            f.write('-' * len('Panel %d:' % i))
            f.write('\n')
            f.write('number of samples: %d\n' % len(means[j]))
            f.write('\n')
            for k in panels[j]:
                f.write(' * %s\n' % k)
            f.write('\n')
            # draw number of events figure
            fig = plt.figure()
            z = np.array(nevents[j])
            ax = fig.add_subplot(1, 1, 1)
            ax.hist(z, bins=10, histtype='step')
            ax.set_xlabel('Number of events')
            plt.tight_layout()
            fig.savefig('nevents_panel_%d.png' % i)
            f.write('![distribution of number of events](nevents_panel_%d.png)\n' % i)
            # draw distribution of means
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            mean = np.array(means[j]).squeeze()
            if len(means[j]) > 1:
                bp = ax.boxplot(np.log10(mean))
                set_color(bp)
            else:
                ax.plot(np.arange(len(mean)), np.log10(mean), 'b+')
                ax.set_xticks(np.arange(len(mean)))
                ax.set_xlim((-1, len(mean)))
            ax.set_xticklabels(panels[j], rotation=90)
	    ax.set_ylabel(r'log_10()')
	    ax.set_title("Distribution of Means")
	    plt.subplots_adjust(bottom=0.8)
	    plt.tight_layout()
            fig.savefig('dist_panel_%d.png' % i)
            f.write('![distribtuion of means by channel](dist_panel_%d.png)\n' % i)
            f.write('\n')
        f.write('Common Markers:\n')
        f.write('---------------\n')
        for j in common:
            f.write(' * %s\n' % j)

        f.write('\n')
        f.write('Laser overview:\n')
        f.write('---------------\n')
        for i in laserabs:
            markers = list(set(laserabs[i]))
            markers.sort()
            f.write('%s\t: %s\n\n' % (i, ', '.join(markers)))
        f.write('\n')
	print 'DONE'


if __name__ == '__main__':
    import os
    os.chdir('/home/jolly/Projects/provide/data/provide_cytof_data_011413')
    report('*/*.fcs', None)
