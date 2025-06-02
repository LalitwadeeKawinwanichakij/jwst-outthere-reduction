# Please see our requirements.txt for the packages required to be installed
#https://drive.google.com/file/d/1PZbyYq-5XxHNNg13UoMmRcjXk3GUeYCx/view?usp=share_link
import argparse
from multiprocessing import Pool
import matplotlib
import pickle
matplotlib.use('Agg')
import astropy.wcs as pywcs
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval
import matplotlib.colors as mcolors
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
from scipy import stats
import glob,os
import warnings
import copy,sys

import matplotlib.image as mpimg
warnings.filterwarnings('ignore')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
myfontsize=22
plt.rcParams.update({'font.size': myfontsize})
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,ScalarFormatter, NullFormatter,MaxNLocator, NullLocator,LogLocator
from grizli.utils import get_line_wavelengths 
from grizli import utils
import matplotlib.patheffects as path_effects
from matplotlib.lines import Line2D
#from matplotlib.patches import Patch
import matplotlib.patches as patches
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d

from scipy.signal import find_peaks
from scipy.ndimage import binary_dilation,center_of_mass
from grizli import multifit
lambda_minus = {'F115W': 1.013, 'F150W': 1.330, 'F200W': 1.751}
lambda_plus = {'F115W': 1.283, 'F150W': 1.671, 'F200W': 2.226}

# Define filter gaps
filter_gaps_NIRISS = [
    (lambda_plus['F115W'], lambda_minus['F150W']),  # gap1
    (lambda_plus['F150W'], lambda_minus['F200W']),  # gap2
    (lambda_plus['F200W'], float('inf'))           # gap3
]



#Taking color scheme from Grizli
grism_colors={'F115W': (0.0, 0.6196078431372549, 0.45098039215686275),
              'F150W': (0.8352941176470589, 0.3686274509803922, 0.0),
              'F200W': (0.8, 0.4745098039215686, 0.6549019607843137)}

parser = argparse.ArgumentParser(description="Make PDF of biopage of each galaxy.")
parser.add_argument("--field", required=True, help="Field name (e.g., sex-00)")
parser.add_argument('--output_dir',type=str,required=True)
parser.add_argument('--cpus',type=int,default=8)
parser.add_argument('--useprior',type=int,default=0)
args = parser.parse_args()
root=args.field
outputdir = args.output_dir
snr_thresh = 7
ew_thresh=50
matched_contam_lambda_thresh=0.02
ncpu = args.cpus
grizli_ewprior_dir='/fred/oz408/NISPureParallel/FIELDS/%s/RedshiftFitting' % root
grizli_extract = '/fred/oz408/NISPureParallel/FIELDS/%s/Extractions' % root
if args.useprior==0:
    grizli_outputdir = grizli_extract
if args.useprior==1:
    grizli_outputdir = grizli_ewprior_dir
cat0=Table.read('%s/%s_phot_apcorr.fits' % (grizli_extract,root)) # all sources 
cat=cat0


#outputdir='/Users/lkawinwanichakij/NISPureParallel/FIELDS/%s' % root
plots_dir= outputdir #'%s/Summary_Plots' % (outputdir)
try:
    os.makedirs(plots_dir, exist_ok=True)
    print(f"Directory '{plots_dir}' is ready.")
except Exception as e:
    print(f"Error creating directory: {e}")

args =np.load('%s/fit_args_wNIRISSphot.npy' % (grizli_extract), allow_pickle=True)[0]

def scale_bar(ax, d, dist=1/0.13, text='1"', color='black', flipped=False, fontsize=15,linewidth=2):
    if flipped:
        p0 = d - d / 15. - dist
        p1 = d / 15.
        ax.plot([p0, p0 + dist], [p1, p1], linewidth=linewidth, color=color)
        ax.text(p0 + dist / 2., p1 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
    else:
        p0 = d / 15.
        ax.plot([p0, p0 + dist], [p0, p0], linewidth=linewidth, color=color)
        ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')

def make_random_cmap(ncolors=256, seed=None):
    """
    Make a matplotlib colormap consisting of (random) muted colors.

    A random colormap is very useful for plotting segmentation images.

    Parameters
    ----------
    ncolors : int, optional
        The number of colors in the colormap.  The default is 256.

    seed : int, optional
        A seed to initialize the `numpy.random.BitGenerator`. If `None`,
        then fresh, unpredictable entropy will be pulled from the OS.
        Separate function calls with the same ``seed`` will generate the
        same colormap.

    Returns
    -------
    cmap : `matplotlib.colors.ListedColormap`
        The matplotlib colormap with random colors in RGBA format.
    """
    from matplotlib import colors

    rng = np.random.default_rng(seed)
    hue = rng.uniform(low=0.0, high=1.0, size=ncolors)
    sat = rng.uniform(low=0.2, high=0.7, size=ncolors)
    val = rng.uniform(low=0.5, high=1.0, size=ncolors)
    hsv = np.dstack((hue, sat, val))
    rgb = np.squeeze(colors.hsv_to_rgb(hsv))

    return colors.ListedColormap(colors.to_rgba_array(rgb))

def get_cmap(segmap, background_color='#000000ff', seed=None):
        
        """
        Define a matplotlib colormap consisting of (random) muted
        colors.

        This is useful for plotting the segmentation array.

        Parameters
        ----------
        

        background_color : Matplotlib color, optional
            The color of the first color in the colormap.
            The color may be specified using any of
            the `Matplotlib color formats
            <https://matplotlib.org/stable/tutorials/colors/colors.html>`_.
            This color will be used as the background color (label = 0)
            when plotting the segmentation image. The default color is
            black with alpha=1.0 ('#000000ff').

        seed : int, optional
            A seed to initialize the `numpy.random.BitGenerator`. If
            `None`, then fresh, unpredictable entropy will be pulled
            from the OS. Separate function calls with the same ``seed``
            will generate the same colormap.

        Returns
        -------
        cmap : `matplotlib.colors.ListedColormap`
            The matplotlib colormap with colors in RGBA format.
        """
        
        
        from matplotlib import colors
       
        #ncolors : int
        #    The number of the colors in the colormap.
        ncolors=np.max(segmap)+1 #len(np.unique(segmap))+1
        cmap = make_random_cmap(ncolors, seed=seed)

        #if background_color is not None:
        cmap.colors[0] = colors.to_rgba(background_color)
        #print(cmap.colors)
        return cmap
    
"""
Function for reading emission/absorption lines, their line fluxes (line flux errors), and equivalent widths
given the threshold in signal-to-noise ratio.
"""
def get_lines(full=None,min_snr=5):
    """
    full: output fits file name from grizli (*full.fits)
    min_snr: minimum signal-to-noise ratio of emission lines to be read from file
    """

       
    
    matching_keys = [key for key in full['COVAR'].header if key.startswith('FLUX_')]
    nlines = len(matching_keys)# full['PRIMARY'].header['NUMLINES']
    line_wave_obs_dict={}
    lines_snr_dict={}
    lines_flux_dict={}
    lines_fluxerr_dict={}
    EW16_dict={}
    EW50_dict={}
    EW84_dict={}
    EW_snr_dict={}
    lines_name=[]
    #lin
    _line_wavelengths, _line_ratios = get_line_wavelengths()
    for i in range(nlines):
        lineflux=full['COVAR'].header['FLUX_%s' % str(i).zfill(3)]
        lineflux_err=full['COVAR'].header['ERR_%s' % str(i).zfill(3)]
        comments=(full['COVAR'].header.comments['FLUX_%s'  % str(i).zfill(3) ]).split()
        if lineflux_err>0:
            _snr_line=lineflux/lineflux_err
        if lineflux_err==0:
            _snr_line=-99
    
        ew16=full['COVAR'].header['EW16_%s' % str(i).zfill(3)]
        ew50=full['COVAR'].header['EW50_%s' % str(i).zfill(3)]
        ew84=full['COVAR'].header['EW84_%s' % str(i).zfill(3)]
        if (_snr_line>min_snr) & (ew50>0):
            line_wave_obs_dict['%s' % comments[0]]=_line_wavelengths['%s' % comments[0]][0]*(1+full['ZFIT_STACK'].header['Z_MAP'])
            lines_snr_dict['%s' % comments[0]]=lineflux/lineflux_err
            lines_flux_dict['%s' % comments[0]]=lineflux
            lines_fluxerr_dict['%s' % comments[0]]=lineflux_err
            # compute 16th, 50th, and 84th percentiles of the distribution of equivalent width  
            EW16_dict['%s' % comments[0]]=full['COVAR'].header['EW16_%s' % str(i).zfill(3)]
            EW50_dict['%s' % comments[0]]=full['COVAR'].header['EW50_%s' % str(i).zfill(3)]
            EW84_dict['%s' % comments[0]]=full['COVAR'].header['EW84_%s' % str(i).zfill(3)]
            EW_snr_dict['%s' % comments[0]]=ew50/((ew84-ew16)/2)
            lines_name.append(comments[0])
    
    lines_prop_dicts={'name':lines_name,
                      'wavelength_obs':line_wave_obs_dict,
                      'flux':lines_flux_dict,
                      'flux_err':lines_fluxerr_dict,
                      'snr_line':lines_snr_dict,
                      'snr_ew':EW_snr_dict,
                      'ew_16':EW16_dict,
                      'ew_50':EW50_dict,
                      'ew_84':EW84_dict}
    return lines_prop_dicts
def count_peaks_pz(full_hdul=None):
    ztab = Table(full_hdul['ZFIT_STACK'].data)
    # Extract data
    zgrid = ztab['zgrid']
    logpdf = np.log10(ztab['pdf'])
    min_height = -4 # Minimum probability density for peaks
    #min_prominence = None #0.005  
    min_prominence = 0.5 # Minimum prominence to avoid small fluctuations
    # Detect peaks with height and prominence filtering
    peaks, properties = find_peaks(logpdf, height=min_height , prominence=min_prominence,distance=10)
    n_peaks_pz = len(peaks)
    return n_peaks_pz,np.array(zgrid[peaks])
def get_best_n_solutions_simple(zgrid=None, chi2=None, n=5, dz=0.05):
    # Find minima indices
    minima_indices, _ = find_peaks(-1 * chi2)
    minima_redshifts = zgrid[minima_indices]
    minima_chi_squared = chi2[minima_indices]
    
    # Sort minima by chi² value
    sorted_indices = np.argsort(minima_chi_squared)
    minima_redshifts = minima_redshifts[sorted_indices]
    minima_chi_squared = minima_chi_squared[sorted_indices]
    
    # Filter solutions to ensure uniqueness within dz
    unique_redshifts = []
    unique_chi_squared = []
    
    for z, chi in zip(minima_redshifts, minima_chi_squared):
        # Check if the redshift is too close to an already selected one
        if not any(abs(z - rz) < dz for rz in unique_redshifts):
            unique_redshifts.append(z)
            unique_chi_squared.append(chi)
            if len(unique_redshifts) == n:  # Stop if we reach the desired number of solutions
                break
    
    # If we still don't have enough solutions, add the next best ones
    for z, chi in zip(minima_redshifts, minima_chi_squared):
        if len(unique_redshifts) == n:
            break
        if z not in unique_redshifts:  # Avoid duplicates
            unique_redshifts.append(z)
            unique_chi_squared.append(chi)

    # Convert to numpy arrays
    unique_redshifts = np.array(unique_redshifts)
    unique_chi_squared = np.array(unique_chi_squared)

    # Compute differences with the lowest chi²
    chi2_differences = unique_chi_squared - unique_chi_squared[0]

    return unique_redshifts, unique_chi_squared, chi2_differences

def get_lines_fromcovar(cov_hdu=None,min_snr=5,redshift=None,observed_frame=False):
    from grizli.utils import get_line_wavelengths 
    """
    full: output fits file name from grizli (*full.fits)
    min_snr: minimum signal-to-noise ratio of emission lines to be read from file
    """
    matching_keys = [key for key in cov_hdu.header if key.startswith('FLUX_')]
    nlines = len(matching_keys)# full['PRIMARY'].header['NUMLINES']
    line_wave_obs_dict={}
    lines_snr_dict={}
    lines_flux_dict={}
    lines_fluxerr_dict={}
    EW16_dict={}
    EW50_dict={}
    EW84_dict={}
    EW_snr_dict={}
    lines_name=[]
    #lin
    _line_wavelengths, _line_ratios = get_line_wavelengths()
    for i in range(nlines):
        lineflux=cov_hdu.header['FLUX_%s' % str(i).zfill(3)]
        lineflux_err=cov_hdu.header['ERR_%s' % str(i).zfill(3)]
        comments=(cov_hdu.header.comments['FLUX_%s'  % str(i).zfill(3) ]).split()
        if lineflux_err>0:
            _snr_line=lineflux/lineflux_err
        if lineflux_err==0:
            _snr_line=-99
     
        ew16=cov_hdu.header['EW16_%s' % str(i).zfill(3)]
        ew50=cov_hdu.header['EW50_%s' % str(i).zfill(3)]
        ew84=cov_hdu.header['EW84_%s' % str(i).zfill(3)]
        if (_snr_line>min_snr) & (ew50>0):
            line_wave_obs_dict['%s' % comments[0]]=_line_wavelengths['%s' % comments[0]][0]*(1+redshift)
            lines_snr_dict['%s' % comments[0]]=lineflux/lineflux_err
            lines_flux_dict['%s' % comments[0]]=lineflux
            lines_fluxerr_dict['%s' % comments[0]]=lineflux_err
            # compute 16th, 50th, and 84th percentiles of the distribution of equivalent width  
            EW16_dict['%s' % comments[0]]=cov_hdu.header['EW16_%s' % str(i).zfill(3)] #* (1 + redshift * observed_frame)
            EW50_dict['%s' % comments[0]]=cov_hdu.header['EW50_%s' % str(i).zfill(3)] #* (1 + redshift * observed_frame)
            EW84_dict['%s' % comments[0]]=cov_hdu.header['EW84_%s' % str(i).zfill(3)] #* (1 + redshift * observed_frame)
            EW_snr_dict['%s' % comments[0]]=ew50/((ew84-ew16)/2)
            lines_name.append(comments[0])
    
    lines_prop_dicts={'name':lines_name,
                      'wavelength_obs':line_wave_obs_dict,
                      'flux':lines_flux_dict,
                      'flux_err':lines_fluxerr_dict,
                      'snr_line':lines_snr_dict,
                      'snr_ew':EW_snr_dict,
                      'ew_16':EW16_dict,
                      'ew_50':EW50_dict,
                      'ew_84':EW84_dict}
    return lines_prop_dicts
#def countlines_highsnr(lines_prop_dicts=None,snr_thresh = 5,ew_thresh=50):
#    #nlines_high_snr = np.zeros(len(grizli_table))

    #for i in range(len(grizli_table)):
    #for i in range(3):
#    if 'snr_ew' in lines_prop_dicts:
#        snr_dict = lines_prop_dicts['snr_ew']
#        ew_dict =  lines_prop_dicts['ew_50']
#        linecontam_dict= lines_prop_dicts['linecontam_flag']
#        
#        nlines_high_snr = sum(1 for snr,ew,licon in zip(snr_dict.values(),ew_dict.values(),linecontam_dict.values() ) if  ((snr > snr_thresh) & (ew>ew_thresh) & (ew<1000) & (licon ==0) ))
#    return nlines_high_snr
def is_in_filter_gap(wavelength, gaps):
    """Check if the wavelength falls within any of the specified filter gaps."""
    for gap_start, gap_end in gaps:
        if gap_start < wavelength < gap_end:
            return True
    return False
def find_band_for_wavelength(wavelength, lambda_minus, lambda_plus):
    """
    Determine which filter band a given emission line wavelength belongs to.

    Parameters:
    - wavelength: float, emission line wavelength in microns.
    - lambda_minus: dict, lower cut-off for each filter.
    - lambda_plus: dict, upper cut-off for each filter.

    Returns:
    - band: str, the band the wavelength belongs to (or None if outside all bands).
    """
    for band in lambda_minus.keys():
        if lambda_minus[band] <= wavelength <= lambda_plus[band]:
            return band  # Return the band that matches
    return None  # If no matching band is found

def countlines_highsnr_v0(lines_prop_dicts=None, snr_thresh=5, ew_thresh=50, filter_gaps=filter_gaps_NIRISS):
    """
    Count the number of reliable emission lines based on SNR, EW, contamination, and filter gaps.

    Parameters:
    - lines_prop_dicts: dict, contains information about emission lines (SNR, EW, contamination, wavelength).
    - snr_thresh: float, threshold for signal-to-noise ratio.
    - ew_thresh: float, threshold for equivalent width.
    - filter_gaps: list of tuples, defines filter gaps as (gap_start, gap_end).

    Returns:
    - nlines_high_snr: int, number of reliable emission lines.
    """
    if 'snr_ew' in lines_prop_dicts:
        snr_dict = lines_prop_dicts['snr_ew']
        ew_dict = lines_prop_dicts['ew_50']
        linecontam_dict = lines_prop_dicts['linecontam_flag']
        wavelength_dict = lines_prop_dicts['wavelength_obs']  # Wavelength in Angstroms
        band_found = find_band_for_wavelength(wavelength_emission, lambda_minus, lambda_plus)
        nlines_high_snr = sum(
            1 for snr, ew, licon, wavelength in zip(
                snr_dict.values(), ew_dict.values(), linecontam_dict.values(), wavelength_dict.values()
            )
            if (
                (snr > snr_thresh)
                and (ew > ew_thresh)
                and (ew < 5000)
                and (licon == 0)
                and not is_in_filter_gap(wavelength / 1e4, filter_gaps)  # Convert to microns and check gaps
            )
        )
    else:
        nlines_high_snr = 0  # Return 0 if 'snr_ew' is not in the input dictionary

    return nlines_high_snr

def compute_medianstack_senscurve(mb=None, selected_band='F200W'):
    """
    Compute the median-stacked sensitivity curve for a given band and return an interpolation function.

    Parameters:
    - mb: Object containing spectral beams.
    - selected_band: str, the filter band to process.

    Returns:
    - sens_interp_func: Interpolation function to get median sensitivity for any wavelength.
    - common_wave_micron: The common wavelength grid in microns.
    - median_sens: The median sensitivity values on the grid.
    """
    pupils = np.array([mb.beams[i].grism.pupil for i in range(len(mb.beams))])

    # Select beams that match the chosen band
    selected_indices = np.where(pupils == selected_band)[0]

    # Store all wavelength and sensitivity curves
    wavelengths = []
    sensitivities = []

    # Iterate through selected beams
    for i in selected_indices:
        wave = mb.beams[i].trace_table['wavelength']
        sens = mb.beams[i].trace_table['sensitivity'] / np.max(mb.beams[i].trace_table['sensitivity'])  # Normalize

        wavelengths.append(wave)
        sensitivities.append(sens)

    # Define a common wavelength grid (choose a fine grid covering all beams)
    common_wave = np.linspace(min(min(w) for w in wavelengths), max(max(w) for w in wavelengths), 1000)

    # Interpolate all sensitivity curves onto the common wavelength grid
    interp_sensitivities = np.array([
        np.interp(common_wave, wavelengths[i], sensitivities[i]) for i in range(len(sensitivities))
    ])

    # Compute the median sensitivity at each wavelength
    median_sens = np.median(interp_sensitivities, axis=0)

    # Convert common wavelength grid to microns
    common_wave_micron = common_wave / 1e4

    # Create an interpolation function
    sens_interp_func = interp1d(common_wave_micron, median_sens, bounds_error=False, fill_value="extrapolate")

    return sens_interp_func, common_wave_micron, median_sens
def countlines_highsnr(lines_prop_dicts=None, snr_thresh=5, ew_thresh=50, 
                       filter_gaps=filter_gaps_NIRISS, mb=None, threshold_sensitivity=0.9):
    """
    Count the number of reliable emission lines based on SNR, EW, contamination, filter gaps, 
    and optionally median sensitivity.

    Parameters:
    - lines_prop_dicts: dict, contains information about emission lines (SNR, EW, contamination, wavelength).
    - snr_thresh: float, threshold for signal-to-noise ratio.
    - ew_thresh: float, threshold for equivalent width.
    - filter_gaps: list of tuples, defines filter gaps as (gap_start, gap_end).
    - mb: Object containing spectral beams. If None, median sensitivity check is skipped.
    - threshold_sensitivity: float, minimum allowed sensitivity for an emission line to be considered valid.

    Returns:
    - nlines_high_snr: int, number of reliable emission lines.
    """
    if 'snr_ew' in lines_prop_dicts:
        #snr_dict = lines_prop_dicts['snr_ew']
        snr_dict = lines_prop_dicts['snr_line'] #?
        ew_dict = lines_prop_dicts['ew_50']
        linecontam_dict = lines_prop_dicts['linecontam_flag']
        wavelength_dict = lines_prop_dicts['wavelength_obs']  # Wavelength in Angstroms

        nlines_high_snr = 0  # Counter for valid emission lines

        for snr, ew, licon, wavelength in zip(
            snr_dict.values(), ew_dict.values(), linecontam_dict.values(), wavelength_dict.values()
        ):
            # Convert wavelength to microns
            wavelength_micron = wavelength / 1e4  

            # Find the corresponding band for the emission line
            band_found = find_band_for_wavelength(wavelength_micron, lambda_minus, lambda_plus)

            if band_found is None:
                continue  # Skip if no matching band is found

            # Default to True if no sensitivity check is applied
            passes_sensitivity_check = True  

            # If `mb` is provided, compute median sensitivity and check threshold
            if mb is not None:
                # Compute median sensitivity curve & interpolation function
                sens_interp_func, _, _ = compute_medianstack_senscurve(mb=mb, selected_band=band_found)

                # Get median sensitivity at emission line wavelength
                median_sensitivity_at_lambda = sens_interp_func(wavelength_micron)

                # Update flag based on sensitivity threshold
                passes_sensitivity_check = (median_sensitivity_at_lambda >= threshold_sensitivity)

            # Apply all selection criteria, including sensitivity check if applicable
            if (
                (snr > snr_thresh)
                and (ew > ew_thresh)
                and (ew < 5000)
                and (licon == 0)
                and not is_in_filter_gap(wavelength_micron, filter_gaps)  # Check filter gap
                and passes_sensitivity_check  # Check sensitivity only if mb is provided
            ):
                nlines_high_snr += 1  # Count valid line

    else:
        nlines_high_snr = 0  # Return 0 if 'snr_ew' is not in the input dictionary

    return nlines_high_snr
def compute_integrated_snr(full_hdul=None, linename='Ha', dilation_iteration=10, dilation_iteration_refined=5, show_plot=True,
                          verbose=False):
    """
    Computes the integrated SNR of an emission line map, with an option to refine the centroid position.

    Parameters:
    - full_hdul: FITS HDU list containing emission line maps.
    - linename: str, name of the emission line (e.g., 'Ha', 'OIII').
    - dilation_iteration: int, number of dilation iterations for initial region.
    - dilation_iteration_refined: int, number of dilation iterations for refined region.
    - show_plot: bool, whether to display the emission line map with contours.

    Returns:
    - integrated_snr: float, integrated signal-to-noise ratio of the emission line region.
    """

    # Check if the requested emission line exists in the FITS header
    if ('LINE', linename) not in full_hdul or ('LINEWHT', linename) not in full_hdul:
        print(f"Error: '{linename}' not found in FITS file.")
        return np.nan  # Return NaN if the line does not exist

    # Extract emission line map and weight map
    emission_line_map = full_hdul['LINE', linename].data
    wht_map = full_hdul['LINEWHT', linename].data  # Extract weight map

    # Convert weight map to error map (handling zeros safely)
    with np.errstate(divide='ignore', invalid='ignore'):
        err_map = np.where(wht_map > 0, 1.0 / np.sqrt(wht_map), np.nan)  # Convert weight to error

    # Initial center pixel
    center_x, center_y = 100, 100

    # Create a mask for the central region (initial seed)
    seed_mask = np.zeros_like(emission_line_map, dtype=bool)
    seed_mask[center_y, center_x] = True  # Mark the center pixel

    # Perform morphological dilation to grow the region
    dilated_mask = binary_dilation(seed_mask, iterations=dilation_iteration)  # Expand region

    # Compute flux-weighted centroid within the dilated region
    y_indices, x_indices = np.where(dilated_mask)  # Get pixel coordinates
    flux_values = emission_line_map[dilated_mask]

    # Avoid division by zero
    if np.sum(flux_values) > 0:
        new_center_x = np.sum(x_indices * flux_values) / np.sum(flux_values)
        new_center_y = np.sum(y_indices * flux_values) / np.sum(flux_values)
    else:
        new_center_x, new_center_y = center_x, center_y  # Keep original if flux is zero

    # Convert to integer for indexing
    new_center_x = int(round(new_center_x))
    new_center_y = int(round(new_center_y))

    #print(f'Refined centroid: ({new_center_x}, {new_center_y})')

    # Extract flux and error values from the refined region
    seed_mask_refined = np.zeros_like(emission_line_map, dtype=bool)
    try:
        seed_mask_refined[new_center_y, new_center_x] = True  # Mark refined center pixel
    except:
        pass
    refined_mask = binary_dilation(seed_mask_refined, iterations=dilation_iteration_refined)  # Expand region again
    flux_values = emission_line_map[refined_mask]
    error_values = err_map[refined_mask]

    # Compute total flux
    total_flux = np.sum(flux_values)

    # Compute total error using quadrature sum (error propagation)
    total_error = np.sqrt(np.nansum(error_values**2))  # Avoids NaNs affecting computation

    # Compute integrated SNR
    integrated_snr = total_flux / total_error if total_error > 0 else np.nan

    # Plot the emission map with both original and refined centers (if enabled)
    if show_plot:
        fig, ax = plt.subplots(figsize=(6, 6))
        zscale = ZScaleInterval()
        vmin, vmax = zscale.get_limits(emission_line_map)

        im = ax.imshow(emission_line_map, origin='lower', cmap='magma', vmin=vmin, vmax=vmax)
        ax.contour(refined_mask, colors='cyan', linewidths=1)  # Overlay refined region

        # Labels and colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Flux (arbitrary units)")
        ax.set_xlabel("X Pixel")
        ax.set_ylabel("Y Pixel")
        ax.set_title(f"Emission Line Map: Refined Flux around ({new_center_x:.1f}, {new_center_y:.1f})", fontsize=12)

        plt.show()

    # Print results
    if verbose:
        print(f"Initial center: ({center_x}, {center_y})")
        print(f"Refined centroid: ({new_center_x:.2f}, {new_center_y:.2f})")
        print(f"Total flux in refined region: {total_flux:.3f}")
        print(f"Total error in refined region: {total_error:.3f}")
        print(f"Integrated SNR: {integrated_snr:.3f}")

    return integrated_snr

def get_sci_hdul_extension(stack_hdul=None):
    filters=[stack_hdul[0].header['GRISM%s' % str(i+1).zfill(3)] for i in range(len(stack_hdul[0].header['GRISM0*']))]
    #print(filters)
    sci_ext_dict={}
    for k in filters:
        sci_ext_dict[k] = [] 
    for i in range(len(stack_hdul)):
        try:
            for f in filters:
                if(('SCI' in stack_hdul[i].header['EXTNAME']) & ('PA' in stack_hdul[i].header)  
                   & (stack_hdul[i].header['GRISM']==f)):
                    sci_ext_dict[f].append(i)     
        except:
            pass
    #print(sci_ext_dict)
    return filters,sci_ext_dict
def compute_2dcontam(stack_hdul=None,sci_ext=None,showplot=True,contam_threshold=0.1):
    import copy
    from scipy.ndimage import label, center_of_mass,measurements
    from scipy import ndimage
    from astropy.visualization import ZScaleInterval
    if showplot:
        fig,axs=plt.subplots(3,1,figsize=(14,12),dpi=150)
    #sci_ext=10
    contam_ext=sci_ext+2
    continu_ext= sci_ext+3
    spec2d_contam =stack_hdul[contam_ext].data
    spec2d_continu=stack_hdul[continu_ext].data
    spec2d_err = 1/np.sqrt(stack_hdul[sci_ext+1].data)
    wmin =stack_hdul[sci_ext].header['WMIN']
    wmax = stack_hdul[sci_ext].header['WMAX']
    N=  np.shape(stack_hdul[sci_ext].data)[1]
    Ny= np.shape(stack_hdul[sci_ext].data)[0]
    # #dilate contamination 3x3 structuring element with connectivity 1, used by default
    struct1 = ndimage.generate_binary_structure(2, 1)
    _contam_mask = np.zeros_like(spec2d_contam)
    #thresh_contam = np.mean(spec2d_contam) #+ 1*np.std(spec2d_contam)
    mean_err_spec2d = np.nanmean(spec2d_err[np.isfinite(spec2d_err)])
    thresh_contam = 1.0*mean_err_spec2d
    #print('thresh_contam ',thresh_contam )
    _contam_mask[spec2d_contam>thresh_contam] = 1
    contam_mask = ndimage.binary_dilation(_contam_mask, structure=struct1,iterations=2)

    # detect continuum
    _modelcontinuum_mask = np.zeros_like(spec2d_continu)
    thresh_continuum=np.mean(spec2d_continu) + 1* np.std(spec2d_continu)
    _modelcontinuum_mask[ spec2d_continu>thresh_continuum]= 1
    modelcontinuum_mask  =ndimage.binary_dilation(_modelcontinuum_mask, structure=struct1,iterations=2)
    overlap =  ((contam_mask ==1) & (modelcontinuum_mask==1))
    # compute fraction of contamination on the observed spectra
    #npix_overlap = len(np.where(overlap==1)[0])
    #npix_spectrum = len(np.where(modelcontinuum_mask==1)[0])
    #contam_frac = npix_overlap /npix_spectrum
    area_contam =  np.sum(overlap)
    area_spectrum =np.sum(modelcontinuum_mask)
    contam_frac = area_contam /area_spectrum
    wavelength_centroids= None
    x_centroids=None
    y_centroids=None
    all_contam_blob_areas=0
    if contam_frac>0 : #contam_threshold: 
        labeled_array, num_features = label(overlap)
        #print('number of centroids of contaminations observed spectrum',num_features)
        if num_features>0:
            colors = plt.cm.viridis(np.linspace(0,1,num_features))
            # Calculate centroids of each labeled component
            #centroids = center_of_mass(overlap, labeled_array, range(1, num_features + 1))
            centroids = center_of_mass(overlap*spec2d_contam, labeled_array, range(1, num_features + 1))
            # Calculate the area (number of pixels) for each blob
            blob_areas = measurements.sum(overlap, labeled_array, index=np.arange(1, num_features + 1))

            #x_centroids= np.zeros(len(centroids))
            #y_centroids= np.zeros(len(centroids))
            x_centroids = np.array([])
            y_centroids = np.array([])
            for i,centroid in enumerate(centroids):
                #print(centroid[1], centroid[0])
                dy=centroid[0]-(0.5*Ny)
                blob_fraction = blob_areas[i]/area_spectrum
                if (np.abs(dy)<7) & (blob_fraction>0.05) :
                    if showplot:
                        axs[0].scatter(centroid[1], centroid[0],  s=100,marker='x',color=colors[i],label=r'$A=%.1f(%.2f),dy=%.1f$' % (blob_areas[i],blob_fraction, dy ))  # Centroids are marked with red plus
                    x_centroids=np.hstack([ x_centroids,centroid[1]])
                    y_centroids=np.hstack([y_centroids, centroid[0]])
                    all_contam_blob_areas=all_contam_blob_areas+blob_areas[i]
            if showplot:        
                axs[0].legend(loc=2,ncol=3,fontsize=myfontsize-10)
        # Calculate the wavelength
        if x_centroids is not None:
            wavelength_centroids = wmin + (x_centroids / (N - 1)) * (wmax - wmin)
        #print('wavelength_centroids',wavelength_centroids)
    #if contam_frac<contam_threshold:
    if contam_frac==0:
        wavelength_centroids=None
        

    if showplot:
        
        zmin,zmax =ZScaleInterval().get_limits(stack_hdul[sci_ext].data)
        axs[0].imshow(stack_hdul[sci_ext].data,origin='lower',aspect='auto',vmin=zmin,vmax=zmax)
        axs[1].imshow(stack_hdul[continu_ext].data,origin='lower',aspect='auto',vmin=zmin,vmax=zmax)
        axs[2].imshow(stack_hdul[contam_ext].data,origin='lower',aspect='auto',vmin=zmin,vmax=zmax)
        axs[2].imshow(np.log10(contam_mask*10),origin='lower',aspect='auto',cmap='Blues',alpha=0.5)
        axs[1].imshow(spec2d_continu,origin='lower',aspect='auto',vmin=zmin,vmax=zmax)
        axs[1].imshow(np.log10(modelcontinuum_mask*10),origin='lower',aspect='auto',cmap='Blues',alpha=0.5)
        
        for i in [0,1,2]:
            axs[i].axhline(y=Ny/2,color='cyan',ls='--')
            if i==0:
                axs[i].imshow(np.log10(overlap),origin='lower',aspect='auto',cmap='Reds_r',alpha=0.5)
        
        axs[0].set_title('sci')
        axs[1].set_title('continuum')
        axs[2].set_title('contamination')
    
    return all_contam_blob_areas/area_spectrum,wavelength_centroids,y_centroids
def run_compute_2dcontam(object_id,showplot=True,grizlioutput_dir=None):
    #dir_name='/Users/lkawinwanichakij/NISPureParallel-main/FIELDS/outthere-hudfn/Grizli_redshiftfit_w3DHSTphot'
    stack_hdul=fits.open('%s/%s_%s.stack.fits' % (grizli_extract,root,str(object_id).zfill(5))) 

    filters,sci_ext_dict=  get_sci_hdul_extension(stack_hdul=stack_hdul)
    
    contam2d_frac_dict={}
    contam2d_lambda_dict= {}

    for k in filters:
        contam2d_frac_dict[k] = np.array([] )
        contam2d_lambda_dict[k] = np.array([] )
        
    for k in sci_ext_dict.keys():
        sci_exts = sci_ext_dict[k]
        for j in range(len(sci_exts)):
            #print(sci_exts[j])
            contam_frac,wavelength_centroids,y_centroids= compute_2dcontam(stack_hdul=stack_hdul,sci_ext=sci_exts[j],showplot=showplot,
                                                                contam_threshold=0.1)
            #print(k,j,wavelength_centroids)
            #print('offset from center of image of centroid',k,y_centroids-0.5*Ny)
            #print(contam_frac,wavelength_centroids)
            #contam2d_frac_dict[k].append(contam_frac)
            if wavelength_centroids is not None:
                contam2d_frac_dict[k]=np.hstack([contam2d_frac_dict[k],contam_frac])
                contam2d_lambda_dict[k] = np.hstack([contam2d_lambda_dict[k],wavelength_centroids])
            #contam2d_lambda_dict[k].append(wavelength_centroids)
    return contam2d_frac_dict,contam2d_lambda_dict
 

def make_color_vivid(color, factor=1.5):
    """
    Increase the saturation and value of a color to make it more vivid.
    'factor' > 1 will make the color more vivid.
    """
    # Convert the color to HSV
    hsv_color = mcolors.rgb_to_hsv(mcolors.to_rgb(color))
    
    # Modify saturation and value to make the color more vivid
    hsv_color[1] = min(1.0, hsv_color[1] * factor)  # Saturation
    hsv_color[2] = min(1.0, hsv_color[2] * factor)  # Value (brightness)
    
    # Convert back to RGB
    vivid_color = mcolors.hsv_to_rgb(hsv_color)
    
    return vivid_color

def complementary_color(normalized_color):
    # Calculate the complementary color by subtracting from 1.0
    return tuple(1.0 - c for c in normalized_color)
#import some useful functions


def get_2dspectra_v2(hdu=None,axs_f115w=None,
              axs_f150w=None,axs_f200w=None,cmap='viridis_r',use_zscale=False,fs_offset=0):
    #fs_offset = 0# subtitle fontsize offset
    h0 = hdu[0].header 
    NX = h0['NGRISM']
    NY = 0
  
    grisms = OrderedDict()
    axs_allbands=[axs_f115w,axs_f150w,axs_f200w]
    for ig in range(NX):
        g = h0['GRISM{0:03d}'.format(ig+1)]
        NY = np.maximum(NY, h0['N'+g])
        grisms[g] = h0['N'+g]

    for g in grisms: # loop over filter
        if g =='F115W':
            ig=0
        if g =='F150W':
            ig=1
        if g=='F200W':
            ig=2
        iters=range(grisms[g])

        sci_i = hdu['SCI', g]
        wht_i = hdu['WHT', g]

        model_i = hdu['MODEL', g]

        #kern_i = hdu['KERNEL', g]
        h_i = sci_i.header
        clip = wht_i.data > 0
        if clip.sum() == 0:
            clip = np.isfinite(wht_i.data)

        if use_zscale is False:
            avg_rms = 1/np.median(np.sqrt(wht_i.data[clip]))
            vmax = np.maximum(1.1*np.nanpercentile(sci_i.data[clip], 98),
                         5*avg_rms)
            vmin=-0.1*vmax
        if use_zscale:
            zscl = ZScaleInterval(contrast=0.25)
            vmin,vmax= zscl.get_limits(sci_i.data[clip])
        # Spectrum
        sh = sci_i.data.shape
        extent = [h_i['WMIN'], h_i['WMAX'], 0, sh[0]]
        #sci 2D spectra (all PAs) (background subtracted by default)
        axs_allbands[ig][0,2].imshow(sci_i.data,origin='lower', interpolation='Nearest',
                  vmin=-0.1*vmax, vmax=vmax, cmap=cmap,
                  extent=extent, aspect='auto')

        axs_allbands[ig][0,2].set_title('%s' % (sci_i.header['EXTVER']),fontsize=myfontsize-fs_offset)
        #model continuum 2D spectra (all PAs)
        axs_allbands[ig][2,2].imshow(model_i.data,origin='lower', interpolation='Nearest',
                  vmin=vmin, vmax=vmax, cmap=cmap,
                  extent=extent, aspect='auto')       
        # continuum subtracted and contamination subtracted
        #axs_allbands[ig][3,2].imshow(sci_i.data-model_i.data,origin='lower', interpolation='Nearest',
        #          vmin=-0.1*vmax, vmax=vmax, cmap=cmap,
        #          extent=extent, aspect='auto')
        axs_allbands[ig][3,2].axis('off')
        for k in [0,1,2,3]:
                axs_allbands[ig][k,2].set_yticklabels([])
                #axs_allbands[ig][k,2].set_xticklabels([])
                axs_allbands[ig][k,2].xaxis.set_major_locator(MultipleLocator(0.1))

        axs_allbands[ig][1,2].axis('off')
        all_PAs=[]
        for col,ip in enumerate(iters): # loop over pa for each filter
            #print(ip, ig)
            pa = h0['{0}{1:02d}'.format(g, ip+1)]
            all_PAs.append(pa)
            #print(ip, g,pa)
            sci_i = hdu['SCI', '{0},{1}'.format(g, pa)]
            wht_i = hdu['WHT', '{0},{1}'.format(g, pa)]
            contam_i = hdu['CONTAM', '{0},{1}'.format(g, pa)]
            model_i = hdu['MODEL', '{0},{1}'.format(g, pa)]
            h_i = sci_i.header
            sh = sci_i.data.shape
            extent = [h_i['WMIN'], h_i['WMAX'], 0, sh[0]]
            #print('pa',pa)

            # sci 2d spectra (not contamination subtract )
            axs_allbands[ig][0,col].imshow(sci_i.data,origin='lower', interpolation='Nearest',
                  vmin=vmin, vmax=vmax, cmap=cmap,
                  extent=extent, aspect='auto')


            axs_allbands[ig][0,col].set_title('%s' % (sci_i.header['EXTVER']),fontsize=myfontsize-fs_offset)
            # model contamination 
            axs_allbands[ig][1,col].imshow(contam_i.data,origin='lower', interpolation='Nearest',
                  vmin=vmin, vmax=vmax, cmap=cmap,
                  extent=extent, aspect='auto')

            #model continuum
            axs_allbands[ig][2,col].imshow(model_i.data,origin='lower', interpolation='Nearest',
                  vmin=vmin, vmax=vmax, cmap=cmap,
                  extent=extent, aspect='auto')

            # sci - contamination - continuum model (left with emission line, if any)
            axs_allbands[ig][3,col].imshow(sci_i.data-contam_i.data-model_i.data,origin='lower', interpolation='Nearest',
                  vmin=vmin, vmax=vmax, cmap=cmap,
                  extent=extent, aspect='auto')            
        axs_allbands[ig][0,0].text(-0.2,0.4,'science',transform=axs_allbands[ig][0,0].transAxes,fontsize=myfontsize-fs_offset,rotation=90) #  
        axs_allbands[ig][1,0].text(-0.2,0.0,'contam',transform=axs_allbands[ig][1,0].transAxes,fontsize=myfontsize-fs_offset,rotation=90) #  
        axs_allbands[ig][2,0].text(-0.2,0.0,'continuum',transform=axs_allbands[ig][2,0].transAxes,fontsize=myfontsize-fs_offset-2,rotation=90) #
        axs_allbands[ig][3,0].text(-0.3,0.,'continuum'+'\n'+'subtracted',transform=axs_allbands[ig][3,0].transAxes,fontsize=myfontsize-fs_offset-4,rotation=90) #

        for kk,pa in enumerate(all_PAs):
            #print('kk,pa',kk,pa)
            #try:
            #_=np.shape(hdu['SCI','%s,%s' % (g,pa)])
            for jj in [0,1,2,3]:
                axs_allbands[ig][jj,kk].set_yticklabels([])
                axs_allbands[ig][jj,kk].xaxis.set_major_locator(MultipleLocator(0.1))
        if len(all_PAs)<2:
            for jj in [0,1,2,3]:
                axs_allbands[ig][jj,1].axis('off')
def plot_summary(object_id):
    #read 1D spectra 
    import numpy as np
    
    spec1d_hdul = fits.open('%s/%s_%s.1D.fits' % (grizli_outputdir,root,str(object_id).zfill(5)))
    #spec1d_hdul.info()
    filters = [spec1d_hdul[0].header['GR*'][i] for i in range(len(spec1d_hdul[0].header['GR*']))]
    #print('Observed filters:', filters)

    # read grizli's spectral fitting output
    full_hdul = fits.open('%s/%s_%s.full.fits' % (grizli_outputdir,root,str(object_id).zfill(5)))
    
    # read stacked 2D spectra
    spec2d_hdul = fits.open('%s/%s_%s.stack.fits' % (grizli_extract,root,str(object_id).zfill(5)))

    mb_file = '%s/%s_%s.beams.fits' % (grizli_extract,root,str(object_id).zfill(5))
    
    #mb = multifit.MultiBeam(
    #        mb_file,
    #        fcontam=args['fcontam'],
    #        group_name=args['group_name'],
    #        MW_EBV=args['MW_EBV'],
    #        sys_err=args['sys_err'],
    #        verbose=False,
    #        psf=args['use_psf'],
    #        min_mask=args['min_mask'],
    #        min_sens=args['min_sens'],
    #        mask_resid=args['mask_resid'],
    #       )

    
    # read cutout image, segmentation map, and individual beams
    beam_hdul = fits.open(mb_file)
    lines_prop_dict = get_lines(full=full_hdul,min_snr=4)
    if  len(lines_prop_dict['name']) >10:
        full_line_list = ['Lya','OII', 'Hb', 'OIII', 'Ha+NII', 'Ha', 'SII', 'SIII', 'MgII']
        filtered_dict = {'name': [line for line in lines_prop_dict['name'] if line in full_line_list],
                             'wavelength_obs': {line: lines_prop_dict['wavelength_obs'][line] for line in lines_prop_dict['name'] if line in full_line_list},
                             'flux': {line: lines_prop_dict['flux'][line] for line in lines_prop_dict['name'] if line in full_line_list},
                             'flux_err': {line: lines_prop_dict['flux_err'][line] for line in lines_prop_dict['name'] if line in full_line_list},
                             'snr_line': {line: lines_prop_dict['snr_line'][line] for line in lines_prop_dict['name'] if line in full_line_list},
                             'snr_ew': {line: lines_prop_dict['snr_ew'][line] for line in lines_prop_dict['name'] if line in full_line_list},
                             'ew_16': {line: lines_prop_dict['ew_16'][line] for line in lines_prop_dict['name'] if line in full_line_list},
                             'ew_50': {line: lines_prop_dict['ew_50'][line] for line in lines_prop_dict['name'] if line in full_line_list},
                             'ew_84': {line: lines_prop_dict['ew_84'][line] for line in lines_prop_dict['name'] if line in full_line_list}}
        lines_prop_dict= filtered_dict
    if len(lines_prop_dict['name']) ==0:
        lines_prop_dict = get_lines(full=full_hdul,min_snr=2)
    
    # compute 2D contamination fraction
 
    contam2d_frac_dict,contam2d_lambda_dict = run_compute_2dcontam(object_id,showplot=False,grizlioutput_dir=grizli_outputdir)
    contam2d_lambda = np.concatenate(list(contam2d_lambda_dict.values()))
    line_wave_micron = np.array([lines_prop_dict['wavelength_obs'][k] for k in lines_prop_dict['snr_line'].keys()])/1e+4
    if len(line_wave_micron)==0:
         lines_prop_dict['linecontam_flag']={}
        
    if len(line_wave_micron)>0:
        #linecontam_flag = np.zeros(len(line_wave_micron)) # flag for each emission line 0 if good, 1 if matched with contamination
        #for name in lines_prop_dict['name']:
        #    lines_prop_dict['linecontam_flag'][name]  = 0
        lines_prop_dict['linecontam_flag']={}
        if len(contam2d_lambda)>0:
           
            for ii in range(len(line_wave_micron)):
                for jj in range(len(contam2d_lambda)):
                    if np.abs( line_wave_micron[ii] - contam2d_lambda[jj])<matched_contam_lambda_thresh:
                        lines_prop_dict['linecontam_flag'][lines_prop_dict['name'][ii]]=1
                        break
                    if np.abs( line_wave_micron[ii] - contam2d_lambda[jj])>matched_contam_lambda_thresh:
                        lines_prop_dict['linecontam_flag'][lines_prop_dict['name'][ii]]=0
        if len(contam2d_lambda)==0:
            for ii in range(len(line_wave_micron)):
                lines_prop_dict['linecontam_flag'][lines_prop_dict['name'][ii]]=0
                        
        #lines_prop_dict['linecontam_flag']=linecontam_flag
    
    nlines_high_snr = countlines_highsnr(lines_prop_dicts=lines_prop_dict,snr_thresh = snr_thresh,ew_thresh=ew_thresh)# ,mb=mb)
   
    
    def get_lineslegend(lines_prop_dicts=None):
        _snr_line= np.array([lines_prop_dicts['snr_line'][k] for k in lines_prop_dicts['snr_line'].keys()])
        _line_name = np.array(lines_prop_dicts['name'])#[lines_prop_dicts['name'][k] for k in lines_prop_dicts['snr_line'].keys()])
        _snr_ew= np.array([lines_prop_dicts['snr_ew'][k] for k in lines_prop_dicts['snr_line'].keys()])
        _line_ew50 = np.array([lines_prop_dicts['ew_50'][k] for k in lines_prop_dicts['snr_line'].keys()])
        _line_ew16 = np.array([lines_prop_dicts['ew_16'][k] for k in lines_prop_dicts['snr_line'].keys()])
        _line_ew84 = np.array([lines_prop_dicts['ew_84'][k] for k in lines_prop_dicts['snr_line'].keys()])
        _line_wave_obs = np.array([lines_prop_dicts['wavelength_obs'][k] for k in lines_prop_dicts['snr_line'].keys()])
        _linecontam_flag =None
        if 'linecontam_flag' in lines_prop_dicts:
            _linecontam_flag = np.array([lines_prop_dicts['linecontam_flag'][k] for k in lines_prop_dicts['snr_line'].keys()])

        linecontam_flag =None
        if len(_snr_line)>0:
            sortinds = np.array(_line_wave_obs).argsort()
            line_name= _line_name[sortinds]
            line_wave_obs =_line_wave_obs[sortinds]
            snr_line = _snr_line[sortinds]
            line_ew50=_line_ew50[sortinds]
            line_ew16=_line_ew16[sortinds]
            line_ew84=_line_ew84[sortinds]
            linecontam_flag=_linecontam_flag[sortinds]
            #selsnrgt5=(np.where(snr_line>5))[0]
            legend_elements=[]
            labelcolor=[]
            for i in range(len(line_name)):
                upper=line_ew84[i]-line_ew50[i]
                lower=line_ew50[i]-line_ew16[i]
                ew_err= (line_ew84[i]-line_ew16[i])/2.0
                snr_ew=line_ew50[i]/ew_err
                #ax.axvline(x=line_wave_obs[i]/10000.0,color='k' )
                #ax.text(line_wave_obs[i]/10000.0,ymin-0.15*ymax,'%s' %line_name[i],fontsize=myfontsize-7,
                #           rotation=90)
                
                if linecontam_flag is None:
                    if (snr_ew> snr_thresh) & (line_ew50[i]>ew_thresh):
                        legend_elements.append(Line2D([0], [0], color='k', lw=1, label=r'\textbf{%s},$S/N_{\mathrm{line}}=%.2f$,$S/N_{EW}=%.2f$,$EW=%.2f_{-%.2f}^{+%.2f}\AA$' %(line_name[i],snr_line[i],snr_ew,line_ew50[i],lower,upper )))
                    else:   
                        legend_elements.append(Line2D([0], [0], color='k', lw=1, label=r'%s,$S/N_{\mathrm{line}}=%.2f$,$S/N_{EW}=%.2f$,$EW=%.2f_{-%.2f}^{+%.2f}\AA$' %(line_name[i],snr_line[i],snr_ew,line_ew50[i],lower,upper )))
                    
                    labelcolor.append('k')
                if linecontam_flag is not None:
                    if linecontam_flag[i]==1:
                        legend_elements.append(Line2D([0], [0], color='lightgrey', lw=1, 
                                                  label=r'%s*,$S/N_{\mathrm{line}}=%.2f$,$S/N_{EW}=%.2f$,$EW=%.2f_{-%.2f}^{+%.2f}\AA$' %(line_name[i],snr_line[i],snr_ew,line_ew50[i],lower,upper )))
                        labelcolor.append('lightgrey')
                    if linecontam_flag[i]==0:
                        legend_elements.append(Line2D([0], [0], color='k', lw=1, 
                                                  label=r'%s,$S/N_{\mathrm{line}}=%.2f$,$S/N_{EW}=%.2f$,$EW=%.2f_{-%.2f}^{+%.2f}\AA$' %(line_name[i],snr_line[i],snr_ew,line_ew50[i],lower,upper )))
                        labelcolor.append('k')
            return legend_elements,labelcolor
        if len(_snr_line)==0:
            return [], [] 
    def plot_1dspec_grism(ax=None,label_each_lines=True,bbox_to_anchor_legend=(0.5, -0.15),legend_fs_offset=8 ):    
        # add label to each emission line?
        if ax is None:
            fig,ax= plt.subplots(figsize=(14,6),dpi=200)
            fig.suptitle(r'ID=%d, $z_{\mathrm{grism}}=%.4f$' % (full_hdul[0].header['ID'],
                                                            full_hdul[0].header['REDSHIFT']),fontsize=myfontsize+5)
        flux_allbands=np.array([])
        wavelength_allbands= np.array([])
        for i in range(len(filters)):
            spec1d = Table(spec1d_hdul[filters[i]].data)
            if 'pscale' in spec1d.colnames:
                pscale= spec1d['pscale']
            if 'pscale' not in spec1d.colnames:
                pscale=1.0
            flux =(spec1d['flux']/spec1d['flat'])/(1.e-19 )/pscale # 1D spectra
            err = (spec1d['err']/spec1d['flat'])/(1.e-19 ) /pscale
            cont=(spec1d['cont']/spec1d['flat'])/(1.e-19 ) # continuum model
            line = (spec1d['line']/spec1d['flat'])/(1.e-19 )
            contam= (spec1d['contam']/spec1d['flat'])/(1.e-19 )/pscale  # contamination model 
            flux_allbands = np.hstack([flux_allbands,flux])
            wavelength_allbands = np.hstack([wavelength_allbands,spec1d['wave']/(1e+4)])
            ax.step(spec1d['wave']/(1e+4),flux,
                    color=grism_colors[filters[i]],where='mid',linestyle='-',label='lines',lw=1,alpha=0.8,zorder=100)

            ax.errorbar(spec1d['wave']/(1e+4),flux,yerr=np.abs(err),
                        ecolor=grism_colors[filters[i]],marker='o',mfc=grism_colors[filters[i]],
                        mec=grism_colors[filters[i]],color='none',ms=0.5,capsize=2,alpha=0.8,zorder=101)
            ax.plot(spec1d['wave']/(1e+4),contam,color=complementary_color(grism_colors[filters[-1]]),ls='--',label='modeled contamination',alpha=1,
                   path_effects=[path_effects.Stroke(linewidth=2, foreground='w'),
                               path_effects.Normal()])

            #color2 = complementary_color(grism_colors[filters[i]])
            color2 =make_color_vivid(grism_colors[filters[i]], factor=2.5)
            #ax.plot(spec1d['wave']/(1e+4),cont,color=color2,ls='-.',label='modeled continuum',
            #       path_effects=[path_effects.Stroke(linewidth=3, foreground='w'),
            #                   path_effects.Normal()])

            ax.plot(spec1d['wave']/(1e+4),line,color=color2,ls='-',lw=2,label='line',path_effects=[path_effects.Stroke(linewidth=3, foreground='w'),
                               path_effects.Normal()],zorder=99,alpha=1)

        xmin = np.nanmin(wavelength_allbands)
        xmax = np.nanmax(wavelength_allbands)
        ymin = np.nanpercentile(flux_allbands,5)
        ymax = np.nanpercentile(flux_allbands,99)
        #if ymax>15:
        #    ymax=15
        #    ymin=-0.05
        ax.set_ylim(ymin,ymax)
        ax.set_xlim(xmin,xmax)
        ax.axhline(y=0,color='k',ls=':')
        ax.set_xlabel(r'Observed wavelength ($\mu$m)')
        ax.set_ylabel(r'$f_{\lambda}[10^{-19}$erg/s/cm$^{-2}/\AA$]')
        ymin,ymax=ax.get_ylim()
        #ax.legend(loc=1,fontsize=myfontsize-5)
        if (label_each_lines==True) & (len(lines_prop_dict['name'])>0):

            for i in range(len(lines_prop_dict['name'])):
                linename= lines_prop_dict['name'][i]
                snr_ew = lines_prop_dict['snr_ew'][lines_prop_dict['name'][i]]
                line_ew50= lines_prop_dict['ew_50'][lines_prop_dict['name'][i]]
                
                if (snr_ew> snr_thresh) & (line_ew50>ew_thresh):
                    lw=2
                    line_text = r'\textbf{%s}'  % lines_prop_dict['name'][i]
                else:
                    lw=0.8
                    line_text = r'%s'  % lines_prop_dict['name'][i]
                    
                if 'linecontam_flag' not in  lines_prop_dict:
                    ax.axvline(x=lines_prop_dict['wavelength_obs'][linename]/10000.0,color='k',lw=lw )
                    ax.text(lines_prop_dict['wavelength_obs'][linename]/10000.0,ymin-0.005*ymax,
                        line_text,fontsize=myfontsize,
                               rotation=90)
                if 'linecontam_flag' in  lines_prop_dict:
                    if  lines_prop_dict['linecontam_flag'][linename]==0:
                        ax.axvline(x=lines_prop_dict['wavelength_obs'][linename]/10000.0,color='k' ,lw=lw)
                        ax.text(lines_prop_dict['wavelength_obs'][linename]/10000.0,ymin-0.005*ymax,
                        line_text,fontsize=myfontsize,
                               rotation=90)
                    if  lines_prop_dict['linecontam_flag'][linename]==1:
                        ax.axvline(x=lines_prop_dict['wavelength_obs'][linename]/10000.0,color='lightgrey' )
                        ax.text(lines_prop_dict['wavelength_obs'][linename]/10000.0,ymin-0.005*ymax,
                        '%s*'  % lines_prop_dict['name'][i],fontsize=myfontsize,
                               rotation=90,color='lightgrey')
                        
                

        """
        Show photometry
        """
        photom_pivot= np.array( [full_hdul[1].header['PHOTL%s' % str(i).zfill(3)] for i in range(len(full_hdul[1].header['PHOTL*']))])

        photom_flam = np.array( [full_hdul[1].header['PHOTF%s' % str(i).zfill(3)] for i in range(len(full_hdul[1].header['PHOTF*']))])
        photom_eflam = np.array( [full_hdul[1].header['PHOTE%s' % str(i).zfill(3)] for i in range(len(full_hdul[1].header['PHOTE*']))])

        if len(photom_pivot)>0:
            selgood= (np.where(photom_eflam>0))[0]
            ax.errorbar(photom_pivot[selgood]/1.e4, photom_flam[selgood]/1.e-19, photom_eflam[selgood]/1.e-19,
                        marker='s', alpha=0.5, color='k', linestyle='None',ms=12)


        legend_elements = [Line2D([0], [0], color='grey', lw=2,ls='-', label='Modeled spectra'),
                           #Line2D([0], [0], color='grey', lw=1,ls=':', label='Modeled continuum'),
                           Line2D([0], [0], color='grey', lw=1,ls='--', label='Modeled contamination')]
        if len(photom_pivot)>0:
            legend_elements.append(Line2D([0], [0], marker='s', color='w', label='NIRISS photometry',
                                  markerfacecolor='k', markersize=15,alpha=0.5))
        leg0= ax.legend(handles=legend_elements, loc=1,fontsize=myfontsize-legend_fs_offset)
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))

        """
        Show legend of emission lines
        """
        #print(lines_prop_dict)
        legend_elements_lines,labelcolor_lines = get_lineslegend(lines_prop_dicts=lines_prop_dict)
        if len(legend_elements_lines)>0:
            ax.legend(handles=legend_elements_lines,loc='upper center', bbox_to_anchor=bbox_to_anchor_legend,
                          fontsize=myfontsize-7,ncol=3,handlelength=0,handletextpad=0,labelcolor=labelcolor_lines)
            ax.add_artist(leg0)
            ax.grid(True, which='both', axis='both', color='gray', linestyle='-', linewidth=1, alpha=0.3)

        return [xmin,xmax],[ymin,ymax]
    def plot_fullsed(ax=None,ylim=None,phot_ms=15):
        if ax is None:
            fig,ax= plt.subplots(figsize=(14,6),dpi=200)
            fig.suptitle(r'ID=%d, $z_{\mathrm{grism}}=%.4f$' % (full_hdul[0].header['ID'],
                                                            full_hdul[0].header['REDSHIFT']),fontsize=myfontsize+5)


        bestfit_templ = Table(full_hdul[3].data)
        ax.plot( (bestfit_templ['wave']*u.Angstrom).to(u.micron).value,bestfit_templ['full']/1.e-19,color='maroon',label='Best-fit SED',
               alpha=0.8)
        ax.set_xlim(0,8)
        """
        Show photometry
        """
        photom_pivot= np.array( [full_hdul[1].header['PHOTL%s' % str(i).zfill(3)] for i in range(len(full_hdul[1].header['PHOTL*']))])

        photom_flam = np.array( [full_hdul[1].header['PHOTF%s' % str(i).zfill(3)] for i in range(len(full_hdul[1].header['PHOTF*']))])
        photom_eflam = np.array( [full_hdul[1].header['PHOTE%s' % str(i).zfill(3)] for i in range(len(full_hdul[1].header['PHOTE*']))])

        if len(photom_pivot)>0:
            selgood= (np.where(photom_eflam>0))[0]
            ax.errorbar(photom_pivot[selgood]/1.e4, photom_flam[selgood]/1.e-19, photom_eflam[selgood]/1.e-19,
                        marker='s', alpha=0.5, color='k', linestyle='None',ms=phot_ms,label='NIRISS photometry')

        #ax.errorbar(photom_pivot[selgood]/1.e4, photom_flam[selgood]/1.e-19, photom_eflam[selgood]/1.e-19,
        #                marker='s', alpha=0.5, color='k', linestyle='None',ms=12,label='3D-HST photometry')
        #ax.set_xscale('log')
        if ylim is None:
            ymin = np.nanpercentile(bestfit_templ['full']/1.e-19,5)
            ymax = np.nanpercentile(bestfit_templ['full']/1.e-19,99)
        if ylim is not None:
            ymin=ylim[0]
            ymax= ylim[1]
        ax.set_ylim(ymin,ymax)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax.legend(loc=1,fontsize=myfontsize-6)
        ax.axhline(y=0,color='k',ls=':')
        ax.set_xlabel(r'Observed wavelength ($\mu$m)')
        ax.set_ylabel(r'$f_{\lambda}[10^{-19}$erg/s/cm$^{-2}/\AA$]')
        ax.grid(True, which='both', axis='both', color='gray', linestyle='-', linewidth=1, alpha=0.3)
        ymin,ymax=ax.get_ylim()
    def plot_pz(axz=None,redshift_hdulist='ZFIT_STACK'):
        redshift_hdulist='ZFIT_STACK'
        if axz is None:
            fig,axz= plt.subplots(figsize=(6,4),dpi=100)
        #try:
        zfit_stack=full_hdul['%s' % redshift_hdulist].data
        Npeaks_pz,z_at_peaks = count_peaks_pz(full_hdul=full_hdul) 
         
        zgrid = zfit_stack['zgrid']
        pdf =zfit_stack['pdf']
        zmi, zma = zgrid.min(), zgrid.max()
        # axz.set_yticks([1,4,9,16,25])
        pzmax = np.log10(pdf.max())
        axz.set_ylim(pzmax-6, pzmax+0.9)
        if (zma-zmi) > 5:
            
            ticks = np.arange(np.ceil(zmi), np.floor(zma), 1)
            lz = np.log(1+zgrid)
            axz.plot(lz, np.log10(pdf), color='maroon', alpha=0.8,zorder=10)
            axz.set_xticks(np.log(1+ticks))
            axz.set_xticklabels(np.cast[int](ticks))
            axz.set_xlim(lz.min(), lz.max())
            #label each z peaks
            if ( Npeaks_pz > 1) & ( Npeaks_pz<5):
                for i in range(len(z_at_peaks)):
                    axz.axvline(np.log(1+z_at_peaks[i]),color='dodgerblue',alpha=0.5,linestyle='--',lw=2,zorder=0)
                    axz.text( np.log(1+z_at_peaks[i]),pzmax+1.0,r'$\mathbf{%.4f}$' % z_at_peaks[i], color='dodgerblue' , rotation =60,
                                  fontsize=myfontsize-4)

                
        else:
            
            axz.plot(zgrid, np.log10(pdf), color='red', alpha=0.5,zorder=10)
            axz.set_xlim(zmi, zma)
            #label each z peaks
            if  Npeaks_pz > 1:
                for i in range(len(z_at_peaks)):
                    axz.axvline(z_at_peaks[i],color='dodgerblue',alpha=0.5,linestyle='--',lw=2,zorder=0)
                    axz.text(z_at_peaks[i],pzmax+1.0,r'$\mathbf{%.4f}$' % z_at_peaks[i], color='dodgerblue' , rotation =90)
       
            
        axz.set_xlabel(r'$z$')
        axz.set_ylabel(r'$\log\ p(z)$'+' / ' + r'$\chi^2=\frac{{{0:.0f}}}{{{1:d}}}={2:.2f}$'.format(full_hdul['%s' % redshift_hdulist].header['CHIMIN'], 
                                                                                                    full_hdul['%s' % redshift_hdulist].header['DOF'],
                                                                                                    full_hdul['%s' % redshift_hdulist].header['CHIMIN']/full_hdul['%s' % redshift_hdulist].header['DOF']),
                      fontsize=myfontsize-5)
       
       
         
        
        axz.grid()
        axz.yaxis.set_major_locator(MultipleLocator(base=1))
    def plot_cutoutdirectim(axs=None,ref_indices=None, ref_filters=None,  ref_filters_list=None,rotate_image=True):
        from scipy.ndimage import rotate
        if axs is None:
            fig,axs= plt.subplots(1,6,figsize=(24,4),dpi=150)
            #fig,axs= plt.subplots(1,len(ref_filters_list)+1,figsize=(4*len(ref_filters_list)+1,4),dpi=150)
        #    axs=axs.ravel()
        #for i in range(len(ref_filters_list)):
        k=0
        for i,filt in enumerate(['F115W', 'F150W', 'F200W']):
            #sel = (np.where(ref_filters==ref_filters_list[i]))[0]
            sel = (np.where(ref_filters==filt))[0]
            if len(sel)>0:
                hdul_ind= ref_indices[sel[0]]
                deltaPix= beam_hdul[hdul_ind].header['PIXSCALE'] # pixel scale
                #print(hdul_ind)
                ref_im=  beam_hdul[hdul_ind].data

                if k==0:
                    z = ZScaleInterval()
                    zmin,zmax = z.get_limits(ref_im)
                #if rotate_image==False:
                pa_v3 = beam_hdul[hdul_ind].header['PA_V3']
                rot_ref_im = rotate(ref_im,pa_v3+90)
                axs[i*2].imshow(rot_ref_im ,origin='lower',vmin=zmin,vmax=zmax) # 
                #if rotate_image==True: # along dispersion direction
                #    pa_v3 = beam_hdul[hdul_ind+2].header['PA_V3']
                #    if 'GR150C' in  beam_hdul[hdul_ind+2].header['FILTER']:
                #        axs[i*2].imshow(rotate(ref_im,pa_v3),origin='lower',vmin=zmin,vmax=zmax)
                #    if 'GR150R' in  beam_hdul[hdul_ind+2].header['FILTER']:
                #        axs[i*2].imshow(rotate(ref_im,pa_v3+90),origin='lower',vmin=zmin,vmax=zmax)

                if beam_hdul[hdul_ind+1].header['PUPIL']=='F200W':
                    seg_im=  beam_hdul[hdul_ind+1].data
                    seg_im_rot=  rotate(seg_im,pa_v3+90,reshape=True, order=0, mode='nearest')
                    seg_im_rot[rot_ref_im==0]=0
                    #seg_im_rot = copy.deepcopy(_seg_im_rot)
                    #seg_im_rot[_seg_im_rot<0] = 0
                    #print(np.min(seg_im),np.max(seg_im),np.min(seg_im_rot),np.max(seg_im_rot))
                    #print(np.unique(seg_im),np.unique(seg_im_rot))
                    segmap_cmap=get_cmap(seg_im_rot)

                    axs[2*i+1].imshow(seg_im_rot,origin='lower',cmap=segmap_cmap)
                    #axs[2*i+1].imshow(seg_im,origin='lower',cmap=get_cmap(seg_im))
                    axs[2*i+1].set_title('Segmentation map(%s)' % (beam_hdul[hdul_ind+1].header['PUPIL']),fontsize=myfontsize-5) 
                    #axs[2*i+1].set_title('Segmentation map',fontsize=myfontsize-8) 
                    axs[2*i+1].axis('off')  
                k=k+1

                axs[2*i].axis('off')
                axs[2*i+1].axis('off')
                axs[2*i].set_title('%s' % filt ,fontsize=myfontsize-5)
                scale_bar(axs[i*2], np.shape(ref_im)[0], dist=1/deltaPix, text=r'$1^{\prime \prime}$',
                      fontsize=myfontsize-2,color='w')
            if len(sel)==0:
                axs[2*i].axis('off')
                axs[2*i+1].axis('off')
    def show_drizzled_lines_v2(axs=None, line_hdu=None,full_line_list=['OII', 'Hb', 'OIII', 'Ha+NII', 'Ha', 'SII', 'SIII'], 
                        size_arcsec=2, cmap='plasma_r',scale=1., dscale=1, 
                        direct_filter=['F140W', 'F160W', 'F125W', 'F105W', 'F110W', 'F098M'],fs_offset=4):
        line_wcs = pywcs.WCS(line_hdu['DSCI'].header)
        pix_size = utils.get_wcs_pscale(line_wcs)
        #pix_size = np.abs(line_hdu['DSCI'].header['CD1_1']*3600)
        majorLocator = MultipleLocator(1.)  # /pix_size)
        N = line_hdu['DSCI'].data.shape[0]/2

        crp = line_hdu['DSCI'].header['CRPIX1'], line_hdu['DSCI'].header['CRPIX2']
        crv = line_hdu['DSCI'].header['CRVAL1'], line_hdu['DSCI'].header['CRVAL2']
        imsize_arcsec = line_hdu['DSCI'].data.shape[0]*pix_size
        # Assume square
        sh = line_hdu['DSCI'].data.shape
        dp = -0.5*pix_size  # FITS reference is center of a pixel, array is edge
        dp = 0
        extent = (-imsize_arcsec/2.-dp, imsize_arcsec/2.-dp,
                  -imsize_arcsec/2.-dp, imsize_arcsec/2.-dp)

        #NL = len(show_lines)

        #xsize = 3*(NL+1)
        #fig = plt.figure(figsize=[xsize, 3.6])

        # Direct
        #ax = fig.add_subplot(1, NL+1, 1)

        dext = 'DSCI'
        # Try preference for direct filter
        for filt in direct_filter:
            if ('DSCI', filt) in line_hdu:
                dext = 'DSCI', filt
                break

        axs[0].imshow(line_hdu[dext].data*dscale, vmin=-0.02, vmax=0.6, cmap=cmap, origin='lower', extent=extent)

        axs[0].set_title('Direct   {0}    z={1:.3f}'.format(line_hdu[0].header['ID'], line_hdu[0].header['REDSHIFT']),
                        fontsize=myfontsize-fs_offset)

        if 'FILTER' in line_hdu[dext].header:
            axs[0].text(0.03, 0.97, line_hdu[dext].header['FILTER'],
                    transform=axs[0].transAxes, ha='left', va='top',fontsize=myfontsize-fs_offset)

        #ax.set_xlabel('RA')
        #ax.set_ylabel('Decl.')

        # Compass
        cosd = np.cos(line_hdu['DSCI'].header['CRVAL2']/180*np.pi)
        dra = np.array([1.5, 1,0,0,0])/3600.*0.12*size_arcsec/cosd
        dde = np.array([0, 0,0,1,1.5])/3600.*0.12*size_arcsec
        cx, cy = line_wcs.all_world2pix(crv[0]+dra, crv[1]+dde, 0)
        cx = (cx-cx.max())*pix_size
        cy = (cy-cy.max())*pix_size
        c0 = 0.95*size_arcsec
        axs[0].plot(cx[1:-1]+c0, cy[1:-1]+c0,
                linewidth=1, color='0.5')
        axs[0].text(cx[0]+c0, cy[0]+c0, r'$E$',
                ha='center', va='center', fontsize=7, color='0.5')
        axs[0].text(cx[4]+c0, cy[4]+c0, r'$N$',
                ha='center', va='center', fontsize=7, color='0.5')

        # 1" ticks
        axs[0].errorbar(-0.5, -0.9*size_arcsec, yerr=0, xerr=0.5, color='k')
        axs[0].text(-0.5, -0.9*size_arcsec, r'$1^{\prime\prime}$', ha='center', va='bottom', color='k')

        # Line maps
        #for i, line in enumerate(show_lines):
            #ax = fig.add_subplot(1, NL+1, 2+i)
        for i in range(len(full_line_list)):
            line = lines_prop_dict['name'][i]
            axs[i+1].imshow(line_hdu['LINE', line].data*scale, vmin=-0.02,
                      vmax=0.6, cmap=cmap, origin='lower', extent=extent)
            axs[i+1].set_title(r'%s %.3f $\mu$m' % (line, line_hdu['LINE', line].header['WAVELEN']/1.e4),
                              fontsize=myfontsize-fs_offset)
            #New compute and show the integrated 2D SNR of each line (Lavender)
            intgrated_snrline = compute_integrated_snr(full_hdul=line_hdu, linename=line,
                                                           dilation_iteration=10, dilation_iteration_refined=5,show_plot=False)
        
            ratio_snr = lines_prop_dict['snr_line'][line]/intgrated_snrline
            color_text='black'
            if ratio_snr< 10 :
                color_text= 'royalblue'
            if ratio_snr>= 10 :
                color_text= 'red'    
            axs[i+1].text(0.1,0.9,r'$\mathbf{SNR_{int,2D} = %.2f}$' % intgrated_snrline,color=color_text,fontsize=myfontsize-fs_offset-4,
                              transform=axs[i+1].transAxes) 

        # End things
        #for ax in fig.axes:
        for i in range(len(axs)):
            axs[i].set_yticklabels([])
            axs[i].set_xticklabels([])
            axs[i].set_xlim(np.array([-1, 1])*size_arcsec)
            axs[i].set_ylim(np.array([-1, 1])*size_arcsec)

            #x0 = np.mean(ax.get_xlim())
            #y0 = np.mean(ax.get_xlim())
            axs[i].scatter(0, 0, marker='+', color='k', zorder=100, alpha=0.5)
            axs[i].scatter(0, 0, marker='+', color='w', zorder=101, alpha=0.5)

            axs[i].xaxis.set_major_locator(majorLocator)
            axs[i].yaxis.set_major_locator(majorLocator)

    
 
    def doplot_z_chi2(extension=None, full_hdul=None, ax=None, target_z=None, chi2_threshold=1):
        extension = 'ZFIT_STACK'  # Must be this HDU so that the full range of zgrid is plotted
        zgrid = full_hdul[extension].data['zgrid']
        chi2 = full_hdul[extension].data['chi2']
        
        try:
            prior = full_hdul[extension].data['prior']
        except:
            prior = None

        pdf = full_hdul[extension].data['pdf']
        zmi, zma = zgrid.min(), zgrid.max()
        ax2 = None
        #print('Hello')
        #print('zmi',zmi)
        #print('zma',zma)
        if (zma - zmi) > 5:
            ticks = np.arange(np.ceil(zmi), np.floor(zma), 1)
            lz = np.log(1 + zgrid)
            # Plot chi² vs. log(1+z) with higher zorder
            ax.plot(lz, chi2, label=r'$\chi^{2}$', color='maroon', alpha=0.7, zorder=3)
            ax.scatter(np.log(1 + full_hdul['ZFIT_BEAM'].header['Z_MAP']),
                       full_hdul['ZFIT_BEAM'].header['CHIMIN'], 
                       marker='x', color='goldenrod', s=100, zorder=4)

            # Plot target redshift line
            if target_z is not None:
                ax.axvline(np.log(1 + target_z), ls='--', alpha=0.5, color='teal', linewidth=2, zorder=2)

            # Find the best redshift solutions
            best_redshifts, best_chi2, chi2_diffs = get_best_n_solutions_simple(zgrid, chi2, n=3)

            # Plot vertical lines for the lowest chi2 solution
            best_z_map = full_hdul['ZFIT_BEAM'].header['Z_MAP']
            ax.axvline(np.log(1 + best_z_map), ls='-', color='black', alpha=0.5, linewidth=2, zorder=1, label='Best Solution')

            # Plot second and third solutions if they meet the chi2 threshold
            for i in range(1, len(best_redshifts)):
                if chi2_diffs[i] < chi2_threshold:
                    ax.axvline(np.log(1 + best_redshifts[i]), ls='--', color='red', alpha=0.5, linewidth=2, zorder=1, label=f'Solution {i+1}')

            # Label the solutions at the top of the plot
            for i, (z_sol, delta_chi2) in enumerate(zip(best_redshifts, chi2_diffs)):
                ax.text(np.log(1 + z_sol), ax.get_ylim()[1] * 0.95, 
                        f"Δχ²={delta_chi2:.2f}\nz={z_sol:.3f}", 
                        color='black' if i == 0 else 'red', 
                        ha='center', fontsize=10, fontweight='bold')

            ax.set_xticks(np.log(1 + ticks))
            ax.set_xticklabels(np.cast[int](ticks))
            ax.set_xlim(lz.min(), lz.max())

            if prior is not None:
                ax2 = ax.twinx()
                ax2.plot(zgrid, np.log(prior), alpha=0.3, color='k')
                ax2.set_ylim(-200, 0)

        else:
            ax.plot(zgrid, chi2, color='maroon', alpha=0.7, zorder=3)
            ax.set_xlim(zmi, zma)
            ax2 = None

        ax.grid(True, alpha=0.4)
        return ax2


    def doplot_z_chi2(extension=None, full_hdul=None, ax=None, target_z=None, chi2_threshold=1):
        extension = 'ZFIT_STACK'  # Must be this HDU so that the full range of zgrid is plotted
        zgrid = full_hdul[extension].data['zgrid']
        chi2 = full_hdul[extension].data['chi2']

        try:
            prior = full_hdul[extension].data['prior']
        except:
            prior = None

        pdf = full_hdul[extension].data['pdf']
        zmi, zma = zgrid.min(), zgrid.max()
        ax2 = None

        if (zma - zmi) > 5:
            ticks = np.arange(np.ceil(zmi), np.floor(zma), 1)
            lz = np.log(1 + zgrid)

            # Plot chi² vs. log(1+z) with higher zorder
            ax.plot(lz, chi2, label=r'$\chi^{2}$', color='maroon', alpha=0.7, zorder=3)
            ax.scatter(np.log(1 + full_hdul['ZFIT_BEAM'].header['Z_MAP']),
                       full_hdul['ZFIT_BEAM'].header['CHIMIN'], 
                       marker='x', color='goldenrod', s=100, zorder=4)

            # Plot target redshift line
            if target_z is not None:
                ax.axvline(np.log(1 + target_z), ls='-', alpha=0.5, color='teal', linewidth=4, zorder=2)

            # Find the best redshift solutions
            best_redshifts, best_chi2, chi2_diffs = get_best_n_solutions_simple(zgrid, chi2, n=3)

            # Plot vertical lines for the lowest chi2 solution
            best_z_map = full_hdul['ZFIT_BEAM'].header['Z_MAP']
            ax.axvline(np.log(1 + best_z_map), ls='-', color='crimson', alpha=0.5, linewidth=2.5, zorder=1, label='Best Solution')

            # Plot second and third solutions if they meet the chi2 threshold
            for i in range(1, len(best_redshifts)):
                if chi2_diffs[i] < chi2_threshold:
                    ax.set_title(r'possibly degenerate with $z=%.2f(\Delta \chi^{2}=%.2f)$' % ( best_redshifts[i],chi2_diffs[i]),
                                     fontsize=10,color='maroon')
                    ax.axvline(np.log(1 + best_redshifts[i]), ls='--', color='lime', alpha=0.5, linewidth=2.5, zorder=1, label=f'Solution {i+1}')
                
            # Label the solutions at the top of the plot
            #for i, (z_sol, delta_chi2) in enumerate(zip(best_redshifts, chi2_diffs)):
            #    ax.text(0,0.95, 
            #            r"$\Delta \chi^2 = {:.2f}$" "\n" r"$z = {:.3f}$".format(delta_chi2, z_sol),  
            #            color='black' if i == 0 else 'red', 
            #            fontsize=10, fontweight='bold',transform=ax.transAxes)

            ax.set_xticks(np.log(1 + ticks))
            ax.set_xticklabels(np.cast[int](ticks))
            ax.set_xlim(lz.min(), lz.max())

            if prior is not None:
                ax2 = ax.twinx()
                ax2.plot(zgrid, np.log(prior), alpha=0.3, color='k')
                ax2.set_ylim(-200, 0)

        else:
            ax.plot(zgrid, chi2, color='maroon', alpha=0.7, zorder=3)
            ax.set_xlim(zmi, zma)
            ax2 = None

        ax.grid(True, alpha=0.4)
        return ax2



    
    def rundoplot_z_chi2(full_hdul=None,ax=None,target_z=None,redshift_hdulist='ZFIT_STACK'):
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        extension='ZFIT_STACK' # Need to be this one to show the entire range of z on x-axis.
        #fig,ax= plt.subplots(figsize=(5,5),dpi=150)
        twin_ax=doplot_z_chi2(extension=redshift_hdulist,full_hdul=full_hdul,ax=ax,target_z=target_z)
       

        fs_offset_inset=10
        fs_offset = 5
        # Add an inset plot
        #inset_ax = inset_axes(ax, width="30%", height="30%", bbox_to_anchor=(0.66,0, 1, 1), bbox_transform=ax.transAxes)  # Adjust size and location
        #_ = doplot_z_chi2(extension='ZFIT_BEAM',full_hdul=full_hdul,ax=inset_ax)
        #for axis,fs in zip([ax,inset_ax], [myfontsize-fs_offset,myfontsize-fs_offset_inset]) :
        for axis,fs in zip([ax], [myfontsize-fs_offset]) :
            axis.set_xlabel(r'$z$',fontsize=fs)   
            axis.set_ylabel(r'$\chi^{2}$',fontsize=fs)   
        ax.tick_params(axis='both', which='major', labelsize=myfontsize-fs_offset)
        if  twin_ax is not None:
            twin_ax.set_ylabel(r'$\ln P(z|F200W)$',fontsize=fs)   
            twin_ax.tick_params(axis='both', which='major', labelsize=myfontsize-fs_offset)
        #inset_ax.tick_params(axis='both', which='major', labelsize=myfontsize-fs_offset_inset)
        #inset_ax.set_title("Zoom-in fit", fontsize=myfontsize-fs_offset_inset)
        #ax.set_yscale('log')
        #ax.plot(zgrid,prior,label=r'$P(z)$')
    def show_spec1d_atz(target_z=None,min_snr=3,ax=None):
        mb_file = '%s/%s_%s.beams.fits' % (grizli_outputdir,root,str(object_id).zfill(5))
        mb = multifit.MultiBeam(mb_file)
        shift, _ = mb.fit_trace_shift( tol=1.0e-3, verbose=False, split_groups=True)
        eazy_phot_file=glob.glob(grizli_outputdir+'/*_eazyphot.pkl')[0] # 
        with open(eazy_phot_file  , "rb") as file:  # Open the file in binary read mode
            int_phot_obj = pickle.load(file)
        phot_i, ii, dd = int_phot_obj.get_phot_dict(mb.ra, mb.dec)
        mb.initialize_masked_arrays()
        if phot_i["flam"] is not None:
            phot = phot_i
        mb.set_photometry(min_err=args['sys_err'], **phot)
        #order=1
        #scale_photometry=order+1
        #try:
        #    scl = mb.scale_to_photometry(
        #        z=0,
        #        method="lm",
        #        order=scale_photometry * 1 - 1,
        #        tol=1.0e-4,
        #        init=None,
        #        fit_background=True,
        #        Rspline=50,
        #        use_fit=True,
        #    )
        #    # tfit=None, tol=1.e-4, order=0, init=None, fit_background=True, Rspline=50, use_fit=True
        #except:
        #    scl = [10.0]
        
        #if hasattr(scl, "status"):
        #    if scl.status > 0:
        #        mb.pscale = scl.x
        #        print('mb.pscale',mb.pscale)
        # use pscale from full.fits
        if 'PSCALEN' in full_hdul[1].header:
            pscale_  =np.zeros(full_hdul[1].header['PSCALEN']+1)
            for i in range(full_hdul[1].header['PSCALEN']+1):   # Hello2
                pscale_[i] =full_hdul[1].header['PSCALE%d' % i]
                mb.pscale = pscale_
        tfit = mb.template_at_z(z=target_z,templates=args['t1'],fit_background=True,fitter=args['fitter'][-1],
                                    bounded_kwargs=args['bounded_kwargs'],use_cached_templates=True)
        oned_hdul = mb.oned_spectrum_to_hdu(tfit=tfit,bin=1,outputfile='%s/%s_%s.1D_test.fits' % (grizli_outputdir,root,str(object_id).zfill(5)),
                                                loglam=args['loglam_1d'],) 
        cov_hdu = fits.ImageHDU(data=tfit["covar"], name="COVAR")
        coeffs_clip = tfit["coeffs"][mb.N :]
        covar_clip = tfit["covar"][mb.N :, mb.N :]
        
        lineEW = utils.compute_equivalent_widths(
                args['t1'], coeffs_clip, covar_clip, max_R=5000, Ndraw=5000, z=tfit["z"]
            )
        #try:
        for ik, key in enumerate(lineEW):
            for j in range(3):
                if not np.isfinite(lineEW[key][j]):
                    lineEW[key][j] = -1.0e30
        #for ik, key in enumerate(tfit["cfit"]):
        #    for j in range(2):
        #        if not np.isfinite( tfit["cfit"][key][j]):
        #             tfit["cfit"][key][j]=-1.0e30

            cov_hdu.header["FLUX_{0:03d}".format(ik)] = tfit["cfit"][key][0], "{0} line flux; erg / (s cm2)".format(key.strip("line "))
            cov_hdu.header["ERR_{0:03d}".format(ik)] = tfit["cfit"][key][1], "{0} line uncertainty; erg / (s cm2)".format(key.strip("line "))

            cov_hdu.header["EW16_{0:03d}".format(ik)] = lineEW[key][0], "Rest-frame {0} EW, 16th percentile; Angstrom".format(key.strip("line "))
            cov_hdu.header["EW50_{0:03d}".format(ik)] = lineEW[key][1], "Rest-frame {0} EW, 50th percentile; Angstrom".format(key.strip("line "))
            cov_hdu.header["EW84_{0:03d}".format(ik)] = lineEW[key][2], "Rest-frame {0} EW, 84th percentile; Angstrom".format(key.strip("line "))
            cov_hdu.header["EWHW_{0:03d}".format(ik)] = (lineEW[key][2] - lineEW[key][0]) / 2, "Rest-frame {0} EW, 1-sigma half-width; Angstrom".format(
                key.strip("line "))
        plot1dspec_from_hdu(spec1d_hdul=oned_hdul,ax=ax,verbose=True)
        #print(cov_hdu.header)
        lines_prop_dict_targetz = get_lines_fromcovar(cov_hdu=cov_hdu,min_snr=2,redshift=target_z)
        ymin,ymax =ax.get_ylim()
        #print('lines_prop_dict_targetz',lines_prop_dict_targetz)
        labeleachline(lines_prop_dict=lines_prop_dict_targetz,ax=ax,ymin=ymin,ymax=ymax)
        #except:
        #    pass
    def plot1dspec_from_hdu(spec1d_hdul=None,ax=None,verbose=True):
        filters = [spec1d_hdul[0].header['GR*'][i] for i in range(len(spec1d_hdul[0].header['GR*']))]
        for i in range(len(filters)):
            spec1d = Table(spec1d_hdul[filters[i]].data)
            if 'pscale' in spec1d.colnames:
                pscale= spec1d['pscale']
            if 'pscale' not in spec1d.colnames:
                pscale=1.0
            #if verbose:
            #    print('pscale',np.array(pscale))
            wave = spec1d['wave']/(1e+4)
            flux =(spec1d['flux']/spec1d['flat'])/(1.e-19 )/pscale # 1D spectra
            err = (spec1d['err']/spec1d['flat'])/(1.e-19 ) /pscale
            cont=(spec1d['cont']/spec1d['flat'])/(1.e-19 ) # continuum model
            line = (spec1d['line']/spec1d['flat'])/(1.e-19 )
            ax.plot(wave,line,color='teal',ls='--',alpha=0.6,lw=4)
    def labeleachline(lines_prop_dict=None,ax=None,ymin=None,ymax=None): 
        for i in range(len(lines_prop_dict['name'])):
            linename= lines_prop_dict['name'][i]
            ax.axvline(x=lines_prop_dict['wavelength_obs'][linename]/10000.0,color='teal',ls=':',alpha=0.6,lw=2 )
            ax.text(lines_prop_dict['wavelength_obs'][linename]/10000.0,ymax+0.005*ymax,
                    '%s'  % lines_prop_dict['name'][i],fontsize=myfontsize,
                           rotation=90,color='teal')

    """
    Make a summary plot
    """
    # Set up a plot 
    # Get indices of FITS hdu which are references images        
    ref_indices =[]
    ref_filters =[]
    for i in range(len(beam_hdul)):
        try:
            if beam_hdul[i].header['EXTNAME']=='REF':
                #print(beam_hdul[i].header['PUPIL'])
                ref_indices.append(i)
                ref_filters.append(beam_hdul[i].header['PUPIL'])
        except:
            pass
    ref_indices= np.array(ref_indices)
    ref_filters = np.array(ref_filters)
    ref_filters_list = np.unique(ref_filters)
    #print('reference bands:',ref_filters_list)


    #fig = plt.figure(constrained_layout=False,dpi=300,figsize=(28,22)) 
    #fig,axs= plt.subplots(1,len(ref_filters_list)+1,figsize=(4*len(ref_filters_list)+1,4),dpi=150)
    fig = plt.figure(constrained_layout=False,dpi=300,figsize=(25,18))
    gs1 = fig.add_gridspec(nrows=11, ncols=14, left=0.05, right=0.95,
                            wspace=0.05,hspace=0.4)
    #gs1 = fig.add_gridspec(nrows=18, ncols=14, left=0.05, right=0.98,
    #                        wspace=0.05,hspace=0.4,bottom=0.01,top=0.92) # Hello
    """
    Direct image cutout
    """
    axs_direct_ims=[]
    #for i in range(len(ref_filters_list)+1):
    for i,filt in enumerate(['F115W', 'F150W', 'F200W']):
        axs_direct_ims.append(fig.add_subplot(gs1[0:2,i*4]))
        axs_direct_ims.append(fig.add_subplot(gs1[0:2,i*4+1]))


    plot_cutoutdirectim(axs=axs_direct_ims,ref_indices=ref_indices, ref_filters=ref_filters,ref_filters_list=ref_filters_list)


    """
    2D spec
    """
    # use four rows
    row=1 #-1
    axs_sci_f115w= [fig.add_subplot(gs1[row+1, 0]),fig.add_subplot(gs1[row+1, 1]),fig.add_subplot(gs1[row+1, 2])]
    axs_sci_f150w= [fig.add_subplot(gs1[row+1,4]),fig.add_subplot(gs1[row+1, 5]),fig.add_subplot(gs1[row+1, 6])]
    axs_sci_f200w= [fig.add_subplot(gs1[row+1,8]),fig.add_subplot(gs1[row+1, 9]),fig.add_subplot(gs1[row+1, 10])]

    axs_cont_f115w= [fig.add_subplot(gs1[row+2, 0]),fig.add_subplot(gs1[row+2, 1]),fig.add_subplot(gs1[row+2, 2])]
    axs_cont_f150w= [fig.add_subplot(gs1[row+2,4]),fig.add_subplot(gs1[row+2, 5]),fig.add_subplot(gs1[row+2, 6])]
    axs_cont_f200w= [fig.add_subplot(gs1[row+2,8]),fig.add_subplot(gs1[row+2, 9]),fig.add_subplot(gs1[row+2, 10])]

    axs_model_f115w= [fig.add_subplot(gs1[row+3, 0]),fig.add_subplot(gs1[row+3, 1]),fig.add_subplot(gs1[row+3, 2])]
    axs_model_f150w= [fig.add_subplot(gs1[row+3,4]),fig.add_subplot(gs1[row+3, 5]),fig.add_subplot(gs1[row+3, 6])]
    axs_model_f200w= [fig.add_subplot(gs1[row+3,8]),fig.add_subplot(gs1[row+3, 9]),fig.add_subplot(gs1[row+3, 10])]
    # # Sci image - continuum
    axs_conti_subtracted_f115w= [fig.add_subplot(gs1[row+4, 0]),fig.add_subplot(gs1[row+4, 1]),fig.add_subplot(gs1[row+4, 2])]
    axs_conti_subtracted_f150w= [fig.add_subplot(gs1[row+4,4]),fig.add_subplot(gs1[row+4, 5]),fig.add_subplot(gs1[row+4, 6])]
    axs_conti_subtracted_f200w= [fig.add_subplot(gs1[row+4,8]),fig.add_subplot(gs1[row+4, 9]),fig.add_subplot(gs1[row+4, 10])]

    # show 2D spectra (science), modeled contamination, and modeled continuum
    axs_f115w_spec2d= np.array([axs_sci_f115w,axs_cont_f115w,axs_model_f115w,axs_conti_subtracted_f115w])
    axs_f150w_spec2d= np.array([axs_sci_f150w,axs_cont_f150w,axs_model_f150w,axs_conti_subtracted_f150w])
    axs_f200w_spec2d= np.array([axs_sci_f200w,axs_cont_f200w,axs_model_f200w,axs_conti_subtracted_f200w])

    get_2dspectra_v2(hdu=spec2d_hdul,axs_f115w=axs_f115w_spec2d,
                  axs_f150w=axs_f150w_spec2d,axs_f200w=axs_f200w_spec2d,fs_offset=8)
    if 'F115W' not in  filters:
        axs=axs_f115w_spec2d.ravel()
        for ax in axs:
            ax.axis('off')
    if 'F150W' not in filters:
        axs=axs_f150w_spec2d.ravel()
        for ax in axs:
            ax.axis('off')
    if 'F200W' not in filters:
        axs=axs_f200w_spec2d.ravel()
        for ax in axs:
            ax.axis('off')
    for axs in [axs_f115w_spec2d,axs_f150w_spec2d,axs_f200w_spec2d]:
        axs=axs.ravel()
        for ax in axs:
            for item in ([ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(fontsize=myfontsize-8)

    """
    Line map
    """
    ncols = len(lines_prop_dict['name'])+1 
    axs_linemap = []
    for i in range(ncols):
        axs_linemap.append(fig.add_subplot(gs1[row+5:row+5+2,i]))
    show_drizzled_lines_v2( axs=axs_linemap,line_hdu=full_hdul, full_line_list=lines_prop_dict['name'],
                          direct_filter=['F200W'],fs_offset=8)   




    phot_ms=8
    bbox_to_anchor_legend=(0.5, -0.18)
    # Create a figure
    #fig = plt.figure(figsize=(18, 4))

    # Define a gridspec with 2 rows and 3 columns
    #gs = fig.add_gridspec(1, 7)

    # Create subplots in different parts of the grid
    row_OT_1d = row+6
    ax1 = fig.add_subplot(gs1[row+7:row+11, 0:2])  
    ax2 = fig.add_subplot(gs1[row+7:row+11, 3:8])
    ax_chi2 = fig.add_subplot(gs1[row+7:row+11, 9:11])
    
    ax3 = fig.add_subplot(gs1[row+7:row+11, 12:])

    

    xlim,ylim = plot_1dspec_grism(ax=ax2,label_each_lines=True,bbox_to_anchor_legend=bbox_to_anchor_legend,
                                 legend_fs_offset=10)
    plot_fullsed(ax=ax1,ylim=None,phot_ms=phot_ms)

    # Define the rectangle parameters
    height = (ylim[1]-ylim[0])  # Height of the rectangle
    x_pos = xlim[0]  # x-coordinate of the lower-left corner
    y_pos = ylim[0]+0.02*height  # y-coordinate of the lower-left corner
    width = xlim[1]-xlim[0]  # Width of the rectangle


    # Create a rectangle patch
    rect = patches.Rectangle((x_pos, y_pos), width, 0.95*height, linewidth=2, ls='--',edgecolor='grey', facecolor='none', alpha=0.8,
                            path_effects=[path_effects.Stroke(linewidth=3, foreground='w'),
                               path_effects.Normal()],zorder=100)

    # Add the rectangle to the plot for showing wavelength range
    ax1.add_patch(rect)
    #fig.subplots_adjust(wspace=0.63)
    #ax2.set_yticklabels([])
    #ax2.set_ylabel('')

    plot_pz(axz=ax3)

    from astropy.coordinates import SkyCoord
    idx = (np.where(cat['id'] == object_id))[0]
    #idx2 = (np.where(grizli_table_w3dhst['id']== object_id))[0]
    #idx3 = (np.where(grizli_table_noprior['id']== object_id))[0]
    try: 
        ra= cat['ra_ot'][idx] #full_hdul[0].header['RA']
    except:
        ra= cat['ra'][idx]
    try:
        dec= cat['dec_ot'][idx]  #full_hdul[0].header['DEC']
    except:
        dec = cat['dec'][idx]
    redchi2 = full_hdul['ZFIT_BEAM'].header['CHIMIN']/full_hdul['ZFIT_BEAM'].header['DOF']
    snr_dict =lines_prop_dict['snr_line']
    selobj2 = (np.where(cat0['id']==object_id))[0]
    try:
        mag_direct = float(-2.5*np.log10(cat0['f200wn_tot_1'][selobj2])+23.9)
        direct_band='F200W'
    except:
        try:
            mag_direct = float(-2.5*np.log10(cat0['f150wn_tot_1'][selobj2])+23.9)
            direct_band='F150W'
        except:
            mag_direct = float(-2.5*np.log10(cat0['f115wn_tot_1'][selobj2])+23.9)
            direct_band='F115W'
            
    # count peak in P(z)
    Npeaks_pz,z_at_peaks = count_peaks_pz(full_hdul=full_hdul)
    npeak_pz_tex=r'$N_{\mathrm{peak,}p(z)}=%d$' % Npeaks_pz 
    # read grizli's output without prior
    if len(glob.glob('/Users/lkawinwanichakij/NISPureParallel-main/FIELDS/outthere-hudfn/Grizli_redshiftfit_no3DHSTphot/%s_%s.full.fits' % (root,str(object_id).zfill(5))) )==1:
        
        full_hdu_noprior= fits.open( '/Users/lkawinwanichakij/NISPureParallel-main/FIELDS/outthere-hudfn/Grizli_redshiftfit_no3DHSTphot/%s_%s.full.fits' % (root,str(object_id).zfill(5)) )
        z_ot_noprior = full_hdu_noprior[1].header['Z_MAP'] #Nancy
        
    else:
        z_ot_noprior =-99
    if len(glob.glob('/Users/lkawinwanichakij/NISPureParallel-main/FIELDS/outthere-hudfn/Grizli_redshiftfit_w3DHSTphot/%s_%s.full.fits' % (root,str(object_id).zfill(5))))==1:
        full_hdu_3dhst = fits.open('/Users/lkawinwanichakij/NISPureParallel-main/FIELDS/outthere-hudfn/Grizli_redshiftfit_w3DHSTphot/%s_%s.full.fits' % (root,str(object_id).zfill(5)))
        z_ot_3dhst = full_hdu_3dhst[1].header['Z_MAP']
    else:
        z_ot_3dhst=-99

    #if (np.abs(cat['z_best'][idx]-full_hdul[0].header['REDSHIFT'])>0.1):
    #show_spec1d_atz(target_z=float(cat['z_dja'][idx]) ,min_snr=3,ax=ax2)
    rundoplot_z_chi2(full_hdul=full_hdul,ax=ax_chi2)#,target_z =float(cat['z_dja'][idx] ))
    
    #print(ra,dec)
    #c= SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    #catalog_3dhst= SkyCoord(ra=cat['ra'],dec=cat['dec'])
    #idx,d2d_3dhst,d3d_3dhst = c.match_to_catalog_3d(catalog_3dhst)
    #z_ot_tex = r'$\mathbf{z_{outthere,grism+phot}=%.4f^{+%.4f}_{-%.4f}}$' % (full_hdul[0].header['REDSHIFT'] ,full_hdul['ZFIT_STACK'].header['Z84']-full_hdul[0].header['REDSHIFT'], full_hdul[0].header['REDSHIFT']- full_hdul['ZFIT_STACK'].header['Z16'])
    #z_ot_tex = r'$\mathbf{z_{outthere,grism+NIRISS~phot,prior}=%.4f^{+%.4f}_{-%.4f}},z_{\mathrm{outthere,grism+NIRISS~phot}}=%.4f, (z_{\mathrm{outthere,grism+3D-HST~phot}}=%.4f)$' % (full_hdul[0].header['REDSHIFT'] ,full_hdul['ZFIT_STACK'].header['Z84']-full_hdul[0].header['REDSHIFT'], full_hdul[0].header['REDSHIFT']- full_hdul['ZFIT_STACK'].header['Z16'],z_ot_noprior,z_ot_3dhst)
    #z_ot_tex = r'$\mathbf{z_{outthere,grism+NIRISS~phot,prior}=%.4f^{+%.4f}_{-%.4f}},z_{\mathrm{outthere,grism+NIRISS~phot}}=%.4f, (z_{\mathrm{outthere,grism+3D-HST~phot}}=%.4f)$' % (full_hdul[0].header['REDSHIFT'] ,full_hdul['ZFIT_BEAM'].header['Z84']-full_hdul[0].header['REDSHIFT'], full_hdul[0].header['REDSHIFT']- full_hdul['ZFIT_BEAM'].header['Z16'],z_ot_noprior,z_ot_3dhst)
    z_ot_tex = r'$\mathbf{z_{outthere,grism+NIRISS~phot}=%.4f^{+%.4f}_{-%.4f}}$' % (full_hdul[0].header['REDSHIFT'] ,full_hdul['ZFIT_BEAM'].header['Z84']-full_hdul[0].header['REDSHIFT'], full_hdul[0].header['REDSHIFT']- full_hdul['ZFIT_BEAM'].header['Z16'])

    # List to hold string components
    string_parts = []
    # Loop through the dictionary to create strings for each key
    for key, values in contam2d_frac_dict.items():
        if len(values) == 0:  # Check if array is empty
            string_parts.append(r'f_{\mathrm{contam,%s}}=0' % key)
        else:
            values_str = ", ".join(r"%.2f"  % v for v in values)  # Format the array values
            string_parts.append(r'f_{\mathrm{contam,%s}}=(%s)' % (key,values_str))
    contam_tex= r"$"+",".join(string_parts)+r"$"
    
    #if d2d_3dhst<0.1*u.arcsec:
        #z_3dhst= cat['z_peak'][idx]
        #qz_3dhst = cat['q_z'][idx]

    #if (cat['z_spec'][idx]>-1):
    #    z_spec_tex=r'$z_{\mathrm{spec}}=%.4f$' % cat['z_spec'][idx]
    #if (cat['z_spec'][idx]==-1):
    #    z_spec_tex=r'$z_{\mathrm{spec}}=-1$'
    #z_best_3dhst_tex= r'$\mathbf{z_{3DHST,best}=%.4f^{+%.4f}_{-%.4f}}$' %  (cat['z_best'][idx],cat['z_best_u68'][idx]-cat['z_best'][idx],cat['z_best'][idx]-cat['z_best_l68'][idx])
    #if cat['z_best_s'][idx] ==1:
    #    z_3dhst_source_tex =r'ground-based spec-z'
    #if cat['z_best_s'][idx] ==2:
    #    z_3dhst_source_tex =r'grism-z'
    #if cat['z_best_s'][idx] ==3:
    #    z_3dhst_source_tex =r'photo-z'
    #if cat['z_best_s'][idx] ==0:
    #    z_3dhst_source_tex =r'star'
    #if (cat['z_max_grism'][idx]>-1):
    #    z_3dhst_max_grism_tex=r'$z_{\mathrm{3D-HST,max~grism}}=%.4f$' % cat['z_max_grism'][idx]
    #if (cat['z_max_grism'][idx]==-1):
    #    z_3dhst_max_grism_tex=r'$z_{\mathrm{3D-HST,max~grism}}=-1$'
    #grism_3dhst_id= str( cat['grism_id'][idx][0]).strip()
    #if  ('G141' in grism_3dhst_id ):
    #    grism_3dhst_id_tex = r'$\mathrm{id}_{\mathrm{3D-HST,grism}}=$%s' % grism_3dhst_id
    #if  ('G141' not in grism_3dhst_id):
    #    grism_3dhst_id_tex =r'$\mathrm{id}_{\mathrm{3D-HST,grism}}=$n/a'
    #threedhst_tex = r'$\mathrm{id}_{\mathrm{3D-HST}}=%d$ %s %s(%s) %s %s' % (cat['id_3dhst'][idx], grism_3dhst_id_tex, z_best_3dhst_tex,z_3dhst_source_tex,
    #                                                           z_3dhst_max_grism_tex, z_spec_tex)
    # add 3dhst 1D spec plot
    #row_3DHST
    #ax_3dhst_linefit = fig.add_subplot(gs1[row+13:,0:6])
    #linefit_png= '../GOODSN_WFC3_V4.1.5_x_outthere0p1arcsec/%s.linefit.png' % (grism_3dhst_id)
    #if len(glob.glob(linefit_png))>0:
    #    image = mpimg.imread(linefit_png)
    #    ax_3dhst_linefit.imshow(image)
    #    ax_3dhst_linefit.axis('off')
    #ax_3dhst_new_zfit = fig.add_subplot(gs1[row+13:,7:]) 
    #new_zfit_png= '../GOODSN_WFC3_V4.1.5_x_outthere0p1arcsec/%s.new_zfit.png' % (grism_3dhst_id)
    #if len(glob.glob(new_zfit_png))>0:
    #    image = mpimg.imread(new_zfit_png)
    #    ax_3dhst_new_zfit.imshow(image)
    #    ax_3dhst_new_zfit.axis('off')

    #dja_tex=r'$z_{\mathrm{spec,DJA}}=%.4f$' % (cat['z_dja'][idx])
    #fig.suptitle(r'object id=%d, ra=%f,dec=%f,$N_{\mathrm{lines~with~SNR}>%d,\mathrm{EW}>%d}=%d, \chi^{2}_{\mathrm{red,outthere}}=%.2f$ %s' % (full_hdul[0].header['ID'],ra,dec,snr_thresh,ew_thresh,nlines_high_snr,redchi2,z_ot_tex) +'\n'+r'mag$_{\mathrm{%s}}=%.1f$, ' % (direct_band,mag_direct)+contam_tex+', '+ nirspec_tex,fontsize=myfontsize-2)
    fig.suptitle(r'object id=%d, ra=%f,dec=%f,$N_{\mathrm{lines~with~SNR}>%d,\mathrm{EW}>%d}=%d, \chi^{2}_{\mathrm{red,outthere}}=%.2f$ %s' % (full_hdul[0].header['ID'],ra,dec,snr_thresh,ew_thresh,nlines_high_snr,redchi2,z_ot_tex)+', '+npeak_pz_tex +'\n'+r'mag$_{\mathrm{%s}}=%.1f$, ' % (direct_band,mag_direct)+contam_tex ,fontsize=myfontsize-2) 

    #if d2d_3dhst>0.1*u.arcsec:  
    #    fig.suptitle(r'object id=%d, ra=%f,dec=%f,$N_{\mathrm{lines~with~SNR}>%d}=%d, \chi^{2}_{\mathrm{red,outthere}}=%.2f$ %s' % (full_hdul[0].header['ID'],ra,dec,snr_thresh,nlines_high_snr,redchi2,z_ot_tex) ,fontsize=myfontsize+2)
    #fig.subplots_adjust(top=0.93,bottom=0.1)
    fig.savefig('%s/%s_%s.pdf' % (plots_dir, root,str(object_id).zfill(5)))
    plt.close()
def run_plot_summary(id):
    
    try:
        
        #filepath = glob.glob(f'{grizli_outputdir}/{root}_{str(id).zfill(5)}.full.fits')
        #local_path = f'{grizli_extract}/{root}_{str(id).zfill(5)}.full.fits'
        #if os.path.isfile(local_path) and not os.path.islink(local_path):
        #print(local_path,os.path.islink(local_path))
        #if len(filepath)>0:
        plot_summary(id)
    except Exception as e:
        print(f"Error processing ID {id}: {e}")
    #finally:
    
if __name__ == '__main__':
    #plot_summary(34) 
    #exit()
    fullfiles = glob.glob('%s/%s_*.full.fits' % (grizli_outputdir, root))
    good_ids=[int(os.path.basename(f).replace(root+'_','').replace('.full.fits','')) for f in fullfiles if not os.path.islink(f)]
    print(len(good_ids))
    if len(glob.glob('%s/%s_*.pdf' % (plots_dir,root)))==0:
        with Pool(processes=ncpu) as pool:
            pool.map(run_plot_summary, good_ids)
            #result = pool.map_async(plot_summary, np.array(cat['id']))
            #_ = result.get()
            #pool.close()
            #pool.join()
        #for i in range(len(cat)):
        #for i in range(len(cat_good3dhstgrism )):
        #    object_id = cat['id'][i]
            #try:
                #if len(glob.glob('%s/%s_%s_result.pdf' % (plots_dir, root,str(object_id).zfill(5))))==0:
        #    try:
        #        plot_summary(object_id =object_id)
        #    except:
        #        pass
    
    if len(glob.glob('%s/%s_*.pdf' % (plots_dir,root)))>0:
        files = glob.glob('%s/%s_*.pdf' % (plots_dir,root))
        done_id = [int(os.path.basename(files[i]).replace('%s_' % root,'').replace('.pdf','')) for i in range(len(files))]
        done_id = np.sort(np.array(done_id))
        all_id = good_ids #np.array(cat['id'])
        notexist_id = list(set(all_id)-set(done_id))
        notexist_id = np.sort(np.array(notexist_id))
        print('number of sources to be run',len(notexist_id))
        print('notexist_id',notexist_id)
        with Pool(processes=ncpu) as pool:
            #result = pool.map_async(plot_summary, notexist_id)
            pool.map(run_plot_summary, notexist_id)
            #_ = result.get()
            #pool.close()
            #pool.join()
        
        #for i in range(len(notexist_id)):
        #for i in [0,1]:
        #    object_id = notexist_id[i]
    
        #    try:
        #        plot_summary(object_id =object_id)
        #    except:
            #    #print(object_id)
        #        pass
                 
             









