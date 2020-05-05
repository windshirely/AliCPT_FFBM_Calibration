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


# Filter requirements.
order = 1
fs = 200       # sample rate, Hz
cutoff = 3.667  # desired cutoff frequency of the filter, Hz
cutoff = 10

# Get the filter coefficients so we can check its frequency response.
b, a = butter_lowpass(cutoff, fs, order)

# Plot the frequency response.
# w, h = freqz(b, a, worN=8000)
# plt.subplot(2, 1, 1)
# plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
# plt.plot(cutoff, 0.5*np.sqrt(2), 'ko')
# plt.axvline(cutoff, color='k')
# plt.xlim(0, 0.5*fs)
# plt.title("Lowpass Filter Frequency Response")
# plt.xlabel('Frequency [Hz]')
# plt.grid()


# Demonstrate the use of the filter.
# First make some data to be filtered.
T = 5.0         # seconds
n = int(T * fs) # total number of samples
t = np.linspace(0, T, n, endpoint=False)
# "Noisy" data.  We want to recover the 1.2 Hz signal from this.
data = np.sin(1.2*2*np.pi*t) + 1.5*np.cos(9*2*np.pi*t) + 0.5*np.sin(12.0*2*np.pi*t)
d0 = np.sin(20*t)
idx = np.arange(10, n, 10)
d0[idx] = 0
data = d0
signal_f = np.fft.rfft(data)
p = np.poly1d([1, 1, 1, 1])
bg = p(t)
bg = bg/np.max(bg)
bg_f = np.fft.rfft(bg)
data = data+bg

# do the fourier transformation
data_f = np.fft.rfft(data)
f = np.fft.rfftfreq(n, d=1. / fs)
# plt.plot(f, np.abs(signal_f), label='signal')
# plt.plot(f, np.abs(bg_f), label='bg')
plt.plot(f, np.abs(data_f), label='data')


# Filter the data, and plot both the original and filtered signals.
y = butter_lowpass_filter(data, cutoff, fs, order)
y_f = np.fft.rfft(y)
plt.plot(f, np.abs(y_f), label='filtered data')
plt.legend()
plt.show()
# plt.subplot(2, 1, 2)
# plt.plot(t, data, 'b-', label='data')
# plt.plot(t, y, 'g-', linewidth=2, label='filtered data')
# plt.xlabel('Time [sec]')
# plt.grid()
# plt.legend()
#
# plt.subplots_adjust(hspace=0.35)
# plt.show()