#!/usr/bin/env python2.7

import numpy
import scipy
import scipy.signal
import pylab
import wave
try:
  import scipy.interpolate
except Exception, e:
  print 'Try installing the latest Numpy-MKL and the latest SciPy.'
  print "If you got those libs from Enthought or your distribution's"
  print "package manager, you may need to find an alternate build or"
  print "else build them from source."
  raise e

import os
import argparse
from math import pi as PI
C = 3e8 # light speed approximation
MOD_PULSE_PERIOD = 20e-3
VCO_FREQ_RANGE = [2400e6, 2591e6] # at 25 degrees, taken from datasheet
#       for my particular VCO given my adjugment of Vtune range.
#
#       MIT's freq range is a default parameter of the RMA
#       if your data filename has the 'mit-' prefix.

# Utility funcs for unit conversion
def meters2feet(meters):
  return meters / 0.3048
def feet2meters(feet):
  return feet * 0.3048

def open_wave(fn):
  '''Returns a tuple of sync_samples, data_samples, and sample_rate.'''
  wavefile = wave.open(fn, 'r')
  raw_data = wavefile.readframes(wavefile.getnframes())
  samples = numpy.fromstring(raw_data, dtype=numpy.int16)
  # try to use 128-bit float for the renormalization
  try:
    normed_samples = samples / numpy.float128(numpy.iinfo(numpy.int16).max)
  except:
    normed_samples = samples / float(numpy.iinfo(numpy.int16).max)
  sync_samples = normed_samples[0::2]
  data_samples = normed_samples[1::2]
  # Need to invert these to match Matlab code for some reason.
  sync_samples *= -1
  data_samples *= -1
  return sync_samples, data_samples, wavefile.getframerate()


# Solely this following function is licensed under the cc-by-sa
# (http://creativecommons.org/licenses/by-sa/3.0/) taken from an answer
# by Joe Kington on Stack Overflow, here:
# http://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array/4495197#4495197
def contiguous_regions(condition):
  """Finds contiguous True regions of the boolean array "condition". Returns
  a 2D array where the first column is the start index of the region and the
  second column is the end index."""
  # Find the indicies of changes in "condition"
  d = numpy.diff(condition)
  idx, = d.nonzero() 
  # We need to start things after the change in "condition". Therefore, 
  # we'll shift the index by 1 to the right.
  idx += 1
  if condition[0]:
    # If the start of condition is True prepend a 0
    idx = numpy.r_[0, idx]
  if condition[-1]:
    # If the end of condition is True, append the length of the array
    idx = numpy.r_[idx, condition.size]
  # Reshape the result into two columns
  idx.shape = (-1,2)
  return idx

def get_sar_frames(sync_samples, data_samples, sample_rate, pulse_period=20e-3):
  '''
  Your sync_samples should look something like this:
  |    ------    ------                   ------    ------
  |----|    |    |    |    |----|----|----|    |    |    |    |----|----
  |         ------    ------                   ------    ------
  Frames are separated by a long enough period of silence between sync oscillations.
  The above simplified ASCII example represents two frames. The length of time between ||s
  is pulse_period (default 20ms), however for the real data there must be at least 250ms of data.
  This restriction is somewhat arbitrary but 150ms is probably the lowest you should go.
  For each frame, this function takes a slice from data_samples at
  the second complete positive-valued period in a given frame, i.e. matching above:
                 |xxxx|
  for the first frame. Then a Hilbert Transform is applied to the frame.
  Finally, the DC phase component is subtracted from all frames,
  and an array of frames is returned.
  '''
  ramp_up_time = pulse_period # the length of the flat top of a sync sample, the time (20 ms) for the frequency modulation to go from lowest to highest.
  minimum_silence_len = sample_rate * 4 * ramp_up_time # arbitrary amount of silence between frames
  # 0.1 is arbitrarily the limit of sensitivity we have for this
  condition = numpy.abs(sync_samples) < 0.1
  silent_regions = contiguous_regions(condition)
  long_enough_silent_regions = filter(lambda x: x[1] - x[0] > minimum_silence_len,
      silent_regions)

  _, start = long_enough_silent_regions.pop(0)
  data_ranges = []
  # we don't need all the samples within a frame, arbitrarily (matching matlab)
  # we just look at the first 0.25 seconds of it. But (different from matlab)
  # because the following code is fragile on the starting points where there
  # may be some extra unwanted oscillations, we set our boundaries a full frame ahead.
  min_time = int(0.25 * sample_rate)
  frame_size = int(0.02 * sample_rate)
  start += frame_size

  for region in long_enough_silent_regions:
    data_ranges.append([start+1, region[0]])
    start = region[1] + frame_size

  sar_frames = []
  for start, end in data_ranges:
    # scan to next positive segment
    if sync_samples[start] > 0: # starting in a positive segment
      break_at_next_positive = False
      for i in xrange(start, min(start + min_time, end)):
        if sync_samples[i] < 0: # hit negative section
          break_at_next_positive = True
        elif sync_samples[i] > 0 and break_at_next_positive:
          start = i
          break
    else: # starting in a negative segment
      for i in xrange(start, min(start + min_time, end)):
        if sync_samples[i] > 0:
          start = i
          break
    # Take Hilbert Transform of each frame.
    # Sadly hilbert() doesn't support float128, so we make a default float
    # type and add the data to it before taking the transform.
    # There's potential for averaging multiple samples per frame instead of just
    # taking one, but that is not done here.
    frame = numpy.zeros(frame_size)
    frame += data_samples[start:start+frame_size]
    frame = scipy.signal.hilbert(frame)
    sar_frames.append(frame)

  # Post-processing: due to the transmit-to-receive antenna coupling, there is
  # a DC phase component in the frames. We can get rid of it by subtracting
  # the mean of all frames from each frame for each frame.
  for i, e in enumerate(sar_frames):
    sar_frames[i] = e - numpy.mean(sar_frames, 0)

  return sar_frames


def RMA(sif, pulse_period=20e-3, freq_range=None, Rs=9.0):
  '''Performs the Range Migration Algorithm.
  Returns a dictionary containing the finished S_image matrix
  and some other intermediary values needed for drawing the image.

  sif is a NxM array where N is the number of SAR frames and M
  is the number of samples within each measurement over the time period
  of frequency modulation increase.

  freq_range should be a tuple of your starting frequency in a range sample and your final frequency.
  If given none, the values from MIT will be used. Please consult your VCO's datasheet data otherwise
  and adjust the constant at the top of this file.

  Rs is distance (in METERS for just this function) to scene center. Default is ~30ft.
  '''
  if freq_range is None:
    freq_range = [2260e6, 2590e6] # Values from MIT

  N, M = len(sif), len(sif[0])

  # construct Kr axis
  delta_x = feet2meters(2/12.0) # Assuming 2 inch antenna spacing between frames.
  bandwidth = freq_range[1] - freq_range[0]
  center_freq = bandwidth/2 + freq_range[0]
  Kr = numpy.linspace(((4*PI/C)*(center_freq - bandwidth/2)), ((4*PI/C)*(center_freq + bandwidth/2)), M)

  # smooth data with hanning window
  sif *= numpy.hanning(M)

  '''STEP 1: Cross-range FFT, turns S(x_n, w(t)) into S(Kx, Kr)'''
  # Add padding if we have less than this number of crossrange samples:
  # (requires numpy 1.7 or above)
  rows = (max(2048, len(sif)) - len(sif)) / 2
  try:
    sif_padded = numpy.pad(sif, [[rows, rows], [0, 0]], 'constant', constant_values=0)
  except Exception, e:
    print "You need to be using numpy 1.7 or higher because of the numpy.pad() function."
    print "If this is a problem, you can try to implement padding yourself. Check the"
    print "README for where to find cansar.py which may help you."
    raise e
  # N may have changed now.
  N = len(sif_padded)

  # construct Kx axis
  Kx = numpy.linspace(-PI/delta_x, PI/delta_x, N)

  freqs = numpy.fft.fft(sif_padded, axis=0) # note fft is along cross-range!
  S = numpy.fft.fftshift(freqs, axes=(0,)) # shifts 0-freq components to center of spectrum

  '''
  STEP 2: Matched filter
  The overlapping range samples provide a curved, parabolic view of an object in the scene. This
  geometry is captured by S(Kx, Kr). Given a range center Rs, the matched filter perfectly
  corrects the range curvature of objects at Rs, partially other objects (under-compsensating
  those close to the range center and overcompensating those far away).
  '''

  Krr, Kxx = numpy.meshgrid(Kr, Kx)
  phi_mf = Rs * numpy.sqrt(Krr**2 - Kxx**2)
  # Remark: it seems that eq 10.8 is actually phi_mf(Kx, Kr) = -Rs*Kr + Rs*sqrt(Kr^2 - Kx^2)
  # Thus the MIT code appears wrong. To conform to the text, uncomment the following line:
  #phi_mf -= Rs * Krr
  # However it is left commented by default because all it seems to do is shift everything downrange
  # closer to the radar by Rs with no noticeable improvement in picture quality. If you do
  # uncomment it, consider just subtracting Krr instead of Krr multiplied with Rs.
  S_mf = S * numpy.exp(1j*phi_mf)

  '''
  STEP 3: Stolt interpolation
  Compensates range curvature of all other scatterers by warping the signal data.
  '''

  kstart, kstop = 73, 108.5 # match MIT's matlab -- why are these values chosen?
  Ky_even = numpy.linspace(kstart, kstop, 1024)

  Ky = numpy.sqrt(Krr**2 - Kxx**2) # same as phi_mf but without the Rs factor.
  try:
    S_st = numpy.zeros((len(Ky), len(Ky_even)), dtype=numpy.complex128)
  except:
    S_st = numpy.zeros((len(Ky), len(Ky_even)), dtype=numpy.complex)
  # if we implement an interpolation-free method of stolt interpolation,
  # we can get rid of this for loop...
  for i in xrange(len(Ky)):
    interp_fn = scipy.interpolate.interp1d(Ky[i], S_mf[i], bounds_error=False, fill_value=0)
    S_st[i] = interp_fn(Ky_even)

  # Apply hanning window again with 1+
  window = 1.0 + numpy.hanning(len(Ky_even))
  S_st *= window

  '''
  STEP 4: Inverse FFT, construct image
  '''

  ifft_len = [len(S_st), len(S_st[0])] # if memory allows, multiply both
  # elements by 4 for perhaps a somewhat better image. Probably only viable on 64-bit Pythons.
  S_img = numpy.fliplr(numpy.rot90(numpy.fft.ifft2(S_st, ifft_len)))

  return {'Py_S_image': S_img, 'S_st_shape': S_st.shape, 'Ky_len': len(Ky), 'delta_x': delta_x, 'kstart': kstart, 'kstop': kstop}

# Based off of example from cansar.py, previously just called out to
# a reduced Matlab/Octave script.
def plot_img(sar_img_data):
  '''Creates the 2D SAR image and saves it as sar_img_data['outfilename'], default sar_image.png.'''
  # Extract S_image, S_st_shape, Ky_len, delta_x, kstart, kstop, Rs, cr1, cr2, dr1, dr2 
  # from sar_img_data
  S_image = sar_img_data['Py_S_image']
  for k, v in sar_img_data.iteritems():
    if k != 'Py_S_image':
      exec('%s=%s' % (k, repr(v)))
  bw = C*(kstop-kstart)/(4*PI)
  max_range = (C*S_st_shape[1]/(2*bw))*1/0.3048

  # data truncation
  dr_index1 = int(round((dr1/max_range)*S_image.shape[0]))
  dr_index2 = int(round((dr2/max_range)*S_image.shape[0]))
  cr_index1 = int(round(S_image.shape[1] * (
    (cr1+Ky_len*delta_x/(2*0.3048)) / (Ky_len*delta_x/0.3048) )))
  cr_index2 = int(round(S_image.shape[1] * (
    (cr2+Ky_len*delta_x/(2*0.3048)) / (Ky_len*delta_x/0.3048) )))

  trunc_image = S_image[dr_index1:dr_index2, cr_index1:cr_index2]
  downrange = numpy.linspace(-1*dr1, -1*dr2, trunc_image.shape[0]) + Rs
  crossrange = numpy.linspace(cr1, cr2, trunc_image.shape[1])

  for i in xrange(0, trunc_image.shape[1]):
    trunc_image[:,i] = (trunc_image[:,i]).transpose() * (abs(downrange*0.3048))**(3/2.0)
  trunc_image = 20 * numpy.log10(abs(trunc_image))

  pylab.figure()
  pylab.pcolormesh(crossrange, downrange, trunc_image, edgecolors='None')
  pylab.plt.gca().invert_yaxis()
  pylab.colorbar()
  pylab.clim([numpy.max(trunc_image)-40, numpy.max(trunc_image)-0])
  pylab.title('Final image')
  pylab.ylabel('Downrange (ft)')
  pylab.xlabel('Crossrange (ft)')
  pylab.axis('equal')
  # Note 'retina' density is about 300, but will increase time of plotting.
  pylab.savefig(sar_img_data['outfilename'], bbox_inches='tight', dpi=200)

def make_sar_image(setup_data):
  '''Gets the frames from an input file, performs the RMA on the SAR data,
  and saves to an output image.'''
  filename = setup_data['filename']

  sif = get_sar_frames(*open_wave(filename), pulse_period=MOD_PULSE_PERIOD)

  if setup_data['bgsub']:
    sif_bg = get_sar_frames(*open_wave(setup_data['bgsub']), pulse_period=MOD_PULSE_PERIOD)
    for i in range(len(sif)):
      if i < len(sif_bg):
        sif[i] -= sif_bg[i]

  Rs = setup_data['Rs']
  freq_range = VCO_FREQ_RANGE
  prefix = filename.split('/')[-1].split('-')[0].lower()
  if prefix == 'mit':
    freq_range = None

  sar_img_data = RMA(sif, pulse_period=MOD_PULSE_PERIOD, freq_range=freq_range, Rs=feet2meters(Rs))

  sar_img_data['outfilename'] = setup_data['outfilename']
  sar_img_data['Rs'] = Rs
  sar_img_data['cr1'] = setup_data['cr1']
  sar_img_data['cr2'] = setup_data['cr2']
  sar_img_data['dr1'] = setup_data['dr1'] + Rs
  sar_img_data['dr2'] = setup_data['dr2'] + Rs

  plot_img(sar_img_data)

def main():
  parser = argparse.ArgumentParser(description="Generate a SAR image outputted by default to 'sar_image.png' from a WAV file of appropriate data.")
  parser.add_argument('-f', nargs='?', type=str, default='mit-towardswarehouse.wav', help="Filename containing SAR data in appropriate format (default: mit-towardswarehouse.wav (prefix filename with 'mit-' to use MIT's frequency range if your VCO range is different))")
  parser.add_argument('-o', nargs='?', type=str, default='sar_image.png', help="Filename to save the SAR image to (default: sar_image.png)")
  parser.add_argument('-rs', nargs='?', type=float, default=30.0, help='Downrange distance (ft) to calibration target at scene center (default: 30)')
  parser.add_argument('-cr1', nargs='?', type=float, default=-80.0, help='Farthest crossrange distance (ft) left of scene center shown in image viewport (default: -80, minimum: -170)')
  parser.add_argument('-cr2', nargs='?', type=float, default=80.0, help='Farthest crossrange distance (ft) right of the scene center shown in image viewport (default: 80, maximum: 170)')
  parser.add_argument('-dr1', nargs='?', type=float, default=1.0, help='Closest downrange distance (ft) away from the radar shown in image viewport (default: 1)')
  parser.add_argument('-dr2', nargs='?', type=float, default=350.0, help='Farthest downrange distance (ft) away from the radar shown in image viewport (default: 350, maximum: 565)')
  parser.add_argument('-bgsub', nargs='?', type=str, default=None, help="Filename containing SAR data representing a background sample that will be subtracted from the main data given by -f (default: None)")

  args = parser.parse_args()

  assert os.path.exists(args.f), "Data file %s not found." % args.f
  try:
    with open(args.o, 'w'):
      pass
  except:
    raise AssertionError('Could not open output file %s for writing.' % args.o)
  assert args.rs > 0, "Rs cannot be 0. It can be 0.0001 or smaller."
  assert (args.cr1 != args.cr2 and -170 <= args.cr1 <= 170 and -170 <= args.cr2 <= 170), "Crossrange values must be between -170 and 170 and not equal."
  assert (args.dr1 != args.dr2 and 1 <= args.dr1 <= 565 and 1 <= args.dr2 <= 565), "Downrange values must be between 1 and 565 and not equal."
  if args.bgsub is not None:
    assert os.path.exists(args.bgsub), "Background substitution file %s not found." % args.bgsub

  setup_data = {'filename': args.f, 'outfilename': args.o, 'Rs': args.rs, 'cr1': args.cr1, 'cr2': args.cr2, 'dr1': args.dr1, 'dr2': args.dr2, 'bgsub': args.bgsub}

  make_sar_image(setup_data)

if __name__ == '__main__':
  main()
