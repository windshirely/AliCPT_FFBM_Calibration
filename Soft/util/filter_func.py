import numpy as np
from scipy.signal import butter, lfilter, freqz
from scipy import signal
import matplotlib.pyplot as plt


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='highpass', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def filter_tod(da, mmp):
    order = 3
    fs = mmp.f_sample
    # cutoff = mmp.srcopt/2
    cutoff = 0.05
    # print(fs, cutoff)
    da = butter_lowpass_filter(da, cutoff, fs, order)
    return da

def filter_tod_v2(da, mmp):
    # win_range = int(mmp.f_sample/(0.5*mmp.srcopt)) # number of sample in one bin
    # win_range = int(mmp.f_sample/mmp.srcopt)+3
    # print(win_range)
    win_range = 8
    win = signal.hann(win_range)
    da = signal.convolve(da, win, mode='same') / sum(win)
    return da


