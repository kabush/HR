## ========================================
## ========================================
##
##  Kayla A. Wilson, BS (2018)
##  Univ. of Arkansas for Medical Sciences
##  Brain Imaging Research Center (BIRC)
##
##  Contribution: Keith A. Bush, Ph.D.
##      - hr trajectory logic mofications
##
## ========================================
## ========================================

import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.io as io
from scipy import signal
from pandas import DataFrame as df
import sys

def find(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]


def calc_hrvs(in_path, input_ids_path, input_times_path, out_path):

    #import kubios file
    beat_times = pd.read_csv(in_path)

    vols = 282
    TR = 2
    Hz = 2000

    # cut kubios file to length of scan
    beat_times = df(beat_times[0:vols * TR * Hz])

    # load in stim times file
    all_stim_times = pd.read_csv(input_times_path, header=None)  # units of seconds

    # load in stim ids file
    stim_ids = pd.read_csv(input_ids_path, header=None)

    # find the ids of image presentations
    ex_ids = np.where(stim_ids == 1)[0]
    next_ids = ex_ids + 1

    # find the times of images presentations
    ex_times = np.array(all_stim_times)[ex_ids[0:45]]

    # find the times of the end of the ITI after image presentation
    next_times = np.array(all_stim_times)[next_ids[0:44]]
    last_time = np.array(ex_times[44] + 8)  # use max possible length for last stimulus
    next_times = np.vstack((next_times, last_time))

    # create a matrix of ranges = [ex_times, next_times] where each row is a range of time
    times = np.hstack((ex_times, next_times))

    # debug...put back in dataframe format to make the syntax work
    times_matrix = df(times, index=list(range(len(times))))

    # for ex_times -1 find beat and beat before that beat, get the mean, & save
    # for ex_times to ex_times + 3 find beats and the beat before ex_times,
    # subtract mu_ibi, & find minimum
    # for next_times - 3 to next_times find beats, subtract mu_ibi, & find maximum
    # save max_dec & min_acc

    # debug ... convert array of arrays to array of values
    beat_times_array = np.squeeze(np.asarray(beat_times.as_matrix()))

    # allocate storage of desired outputs
    max_dec = []
    max_dec_offset = []

    # ----------------------------------------
    # construct 1/2 s trajectories of ibis
    t_intrvs = np.linspace(0.5, 4.0, num=8);
    intrv_hlf = 0.5

    # allocate storage for trajectories
    trajs = np.zeros((len(times_matrix), len(t_intrvs)))

    # construct trajectories
    for i in range(len(times_matrix)):

        ## allocate storage of intermediate outputs
        before_beats = []
        after_beats = []
        after_offsets = []

        ## get current range to search for beats
        row = times_matrix.iloc[i]

        # search for beats occurring within 1.5 sec before ex stimulus
        pre_int = 1.0
        for time in beat_times_array:
            if time > (row.iloc[0] - pre_int) and time < row.iloc[0]:
                before_beats.append(time)

        ## DEBUG: This code fixes bug as this before beat must exist
        ## Therefore, we make worst case assumption about the occurence
        ## of the beat to make before_betas list not empty.
        if not before_beats:
            for time in beat_times_array:
                if time > (row.iloc[0] - 2*pre_int) and time < row.iloc[0]:
                    before_beats.append(time)

        # add additional pre-interval beat to allow for ITI calcs in range
        id_pre = find(beat_times_array, lambda x: x == before_beats[0])

        print(row.iloc[0], row.iloc[1])
        print(before_beats)
        print(id_pre)

        pretime = beat_times_array[id_pre[0] - 1]
        before_beats.insert(0, pretime)

        # search for beats occuring with 3 sec after ex stimulus
        post_int = 6.0
        for time in beat_times_array:
            if time > row.iloc[0] and time < (row.iloc[0] + post_int):
                after_beats.append(time)
                after_offsets.append(time - row.iloc[0])

        # add addtional pre-stimulus beat to all for ITI calcs in range
        id_ex = find(beat_times_array, lambda x: x == after_beats[0])
        pretime_ex = beat_times_array[id_ex[0] - 1]
        after_beats.insert(0, pretime_ex)

        # mean ibi of pre-stimulus beats
        mu_ibi = np.mean(np.diff(before_beats))
        mu_bps = 1/mu_ibi

        # ibi of post-stimulus beats
        dec_ibi = np.diff(after_beats)

        # ----------------------------------------
        # compute HR trajectory post-stimulus
        intrv_bps = []

        print('************************************************')
        print('*** STARTING STIM ***')
        print(dec_ibi)
        print(after_offsets)

        for j in range(len(t_intrvs)):

            intrv = t_intrvs[j]

            intrv_min = intrv - intrv_hlf
            intrv_max = intrv + intrv_hlf
            intrv_width = 2 * intrv_hlf

            ibi_grp = []
            off_grp = []

            for k in range(len(after_offsets)):

                off = after_offsets[k]
                ibi = dec_ibi[k]

                # check for ibi within range and grab if so
                if (off > intrv_min) and (off < intrv_max):
                    off_grp.append(off)
                    ibi_grp.append(ibi)

            # grab one more past end of range to cover end of range
            flag = 0
            for k in range(len(after_offsets)):

                off = after_offsets[k]
                ibi = dec_ibi[k]

                # check for ibi within range and grab if so
                if (off > intrv_max) and (flag == 0):
                    off_grp.append(off)
                    ibi_grp.append(ibi)
                    flag = flag + 1

            print('**************************')
            print(intrv)
            print(off_grp)
            print(ibi_grp)

            # compute the wtd avg of beats-per-second in range of this interval
            # requires coverting ibis (in seconds) to beats-per-second (bps).
            wtd_bps = 0
            for k in range(len(ibi_grp)):

                ibi = ibi_grp[k]
                bps = 1/ibi #convert ibi to bps
                off = off_grp[k]

                print('------------')
                print(off)

                ## if the ibi is first in range then interval
                ## is start of range up to this offset
                if k == 0:
                    print('first element')
                    frac = (off - intrv_min) / intrv_width

                else:

                    # if the beat is past end of range then interval
                    # is only last beat up to end of range
                    if off > intrv_max:
                        print('past max intrv')
                        frac = (intrv_max - off_grp[k - 1]) / intrv_width

                    # if the beat is in range then interval is
                    # the start of the range upto that beat
                    else:
                        print('within intrv not first')
                        frac = (off - off_grp[k - 1]) / intrv_width

                print(frac)
                print(ibi)
                wtd_bps = wtd_bps + frac * bps;

            intrv_bps.append(wtd_bps)

        print(intrv_bps)

        # scale this trajectory of ibis
        intrv_bps_norm = intrv_bps - mu_bps

        # save
        trajs[i, :] = np.asarray(intrv_bps_norm)

    # Write out data
    np.savetxt(out_path + '_t_intrvs.txt', t_intrvs)
    np.savetxt(out_path + '_trajs.txt', trajs)

# ----------------------------------------
# ----------------------------------------

# Gather command line arguments
in_path = sys.argv[1]
input_ids_path = sys.argv[2]
input_times_path = sys.argv[3]
out_path = sys.argv[4]

# Execute
calc_hrvs(in_path, input_ids_path, input_times_path, out_path)
