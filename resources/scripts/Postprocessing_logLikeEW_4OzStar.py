#! /usr/bin/env python
"""
Selected best N redshifts with the lowest chi2, check for the log-likelihood of observing EW at 
that redshift and magnitude of the galaxy. The code will repeat this step until find the solution
with maximum three iterations. If there is no better redshift solution that the intial redshift
from Grizli, the z_best_new will return -99
Note: to make the process faster, we focus on refining redshifts for objects that the first emission line with the highest 
SNR ratio (exceeding 5) are not one of the main emission lines, e.g., main_line_list = ['Ha','Hb','SII','OIII','OII']

The script requires mock spectra with EW are calculated, e.g., 
data_file = '/Users/lkawinwanichakij/JWST/JAGUAR_Mocks/SF_mock_%s_EW.npz' % target_line

"""

import os, glob,pickle,sys,copy,shutil
import scipy
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=False)
#plt.rc('font', family='serif')
#myfontsize=22
#plt.rcParams.update({'font.size': myfontsize})
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,ScalarFormatter, NullFormatter,MaxNLocator, NullLocator,LogLocator
import matplotlib
matplotlib.use('Agg')
import argparse
from astropy.table import Table,vstack
import astropy.io.fits as pyfits
import astropy.units as u
from scipy.stats import gaussian_kde
from scipy.stats import mode

from grizli import utils,prep,multifit,fitting
from grizli.pipeline import auto_script,photoz
import warnings
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser()
parser.add_argument('--field', type=str,required=True)
parser.add_argument('--output_dir',type=str,required=True)

args = parser.parse_args()
root = args.field
new_redshift_dir=args.output_dir
prep_dir= '/fred/oz041/lkawin/NISPureParallel/FIELDS/%s/Prep'  % root # modify if needed
extract_dir= '/fred/oz041/lkawin/NISPureParallel/FIELDS/%s/Extractions'  % root
# new directory to keep new Grizli output.
#the content is the symlinks of initial Grizli outputs for each galaxy)
# the new Grizli outputs will rep# new directory to lace (if the new redshift is found), the existing (symlinks) outputs.
#new_redshift_dir=  '/fred/oz041/lkawin/NISPureParallel/FIELDS/%s/RedshiftFitting'  % root
#os.makedirs(new_redshift_dir, exist_ok=True)
#if len(glob.glob('%s/%s-extracted.fits' % (new_redshift_dir,root)))==0:
#    os.symlink('%s/%s-extracted.fits' % (extract_dir,root),
#            '%s/%s-extracted.fits' % (new_redshift_dir,root))
main_line_list = ['Ha','Hb','SII','OIII','OII']    

"""
Load photometric catalog and initial grizli output
"""
cat=Table.read('%s/%s_phot_apcorr.fits' % (extract_dir,root) )# Aperture corrected photometric catalog of all sources 
grizli_cat=Table.read(extract_dir+'/'+root+'_fitresults.fits') # Grizli initial output catalog.
"""
Load neccessary files for run Grizli
"""

fit_args_file=extract_dir+'/fit_args_wNIRISSphot.npy'
#with open(extract_dir+'/outthere-niriss-%s_eazyphot.pkl' % root, "rb") as f:
#    int_phot_obj = pickle.load(f)
args=  np.load(fit_args_file, allow_pickle=True)[0]
#os.chdir('/fred/oz041/lkawin/NISPureParallel/FIELDS/%s' % root)

"""
Functions to compute line EWs and SNR given the "covar" extension of ***.full.fits
"""
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

"""
Function to compute log likelihood of observed emission line x with EW y at a given redshift and magnitude (F150W or F200W)
"""

def log1p_to_z(log1p_z):
    return np.expm1(log1p_z)  # Inverse of log(1+z)

def z_to_log1p(z):
    return np.log1p(z)  # log(1+z)

def compute_and_plot_joint_density(fig=None,ax=None,target_line=None, specific_ew=None,band=None, 
                                   specific_redshift=None, specific_mag=None,showplot=False):
    """
    Compute and plot the log joint density of EW and log(1+z) for a specific magnitude.

    Parameters:
        target_line (str): Name of the emission line (e.g., "Ha").
        specific_ew (float): Specific EW value for evaluation.
        specific_redshift (float): Specific redshift value for evaluation.
        specific_mag (float): Specific magnitude value for evaluation.
    """
    # Dynamically construct the data file path
    #data_file = f"SF_mock_{target_line}_EW.npz"
    
    data_file = '/fred/oz041/lkawin/NISPureParallel/JAGUAR_Mocks/SF_mock_%s_EW.npz' % target_line
    # Load the mock data
    npz_ew = np.load(data_file)
    #print(npz_ew.files)
    redshifts = npz_ew["redshift"]
    ews = npz_ew["ew"]
    magnitudes = npz_ew["%s_mag" % band]

    # Transform redshifts to log(1+z)
    log1pz = np.log1p(redshifts)

    # Combine data for KDE
    data_points = np.vstack((ews, magnitudes, log1pz))
    kde = gaussian_kde(data_points)

    # Define ranges for EW and log(1+z)
    ew_range = np.linspace(ews.min(), ews.max(), 100)
    log1pz_range = np.linspace(log1pz.min(), log1pz.max(), 100)

    # Generate a grid for evaluation
    ew_grid, log1pz_grid = np.meshgrid(ew_range, log1pz_range)
    grid_points = np.vstack([ew_grid.ravel(), np.full_like(ew_grid.ravel(), specific_mag), log1pz_grid.ravel()])

    # Evaluate the KDE joint density
    joint_density = kde(grid_points).reshape(ew_grid.shape)

    # Convert joint density to log scale (add a small constant to avoid log(0))
    log_joint_density = np.log10(joint_density + 1e-10)

    # Compute log joint density for the specific values
    specific_log1pz = np.log1p(specific_redshift)
    specific_density_point = kde(np.array([[specific_ew], [specific_mag], [specific_log1pz]]))
    specific_log_density = np.log10(specific_density_point[0] + 1e-10)

    # Output the result
    #print(f"Log10 Joint Density at EW={specific_ew},z={log1p_to_z(specific_log1pz):.2f}, Mag={specific_mag:.2f}: {specific_log_density:.6f}")
    print('Log 10 joint density at %s EW = %.1f, z =%.2f, Mag=%.2f: %.3f' % (target_line,specific_ew,specific_redshift,
                                                                            specific_mag,specific_log_density))
    # Plot the joint density as a contour plot
    if showplot:
        if (fig is None) | (ax is None):
            fig, ax = plt.subplots(figsize=(10, 6))
        contour = ax.contourf(ew_grid, log1pz_grid, log_joint_density, levels=30, cmap='viridis')
        cbar = fig.colorbar(contour, ax=ax, label="Log10(Joint Density)")
    
        # Add specific points to the plot
        ax.errorbar(specific_ew, specific_log1pz,color='none', mfc='none',mec='red',marker='*',
                    ms=15, label='%s EW= %.1f at z=%.2f'  % (target_line,specific_ew,log1p_to_z(specific_log1pz) ) )#label=f"EW={specific_ew}, log(1+z)={specific_log1pz:.2f}")
    
        # Customizations
        ax.set_xscale("log")
        ax.set_xlabel("Observed-frame Equivalent Width (EW)")
        ax.set_ylabel("log(1+z)")
        ax.set_title(r'Log Joint Probabilty Density of %s EW and given z and F200W Mag=%.1f' % (target_line,specific_mag))
        ax.legend()
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    return specific_log_density

def get_n_colors(n, cmap_name='viridis'):
    """
    Get a list of N colors from a specified Matplotlib colormap.
    
    Parameters:
    - n (int): Number of colors to generate.
    - cmap_name (str): Name of the Matplotlib colormap (default: 'viridis').
    
    Returns:
    - List of N RGBA color tuples.
    """
    try:
        # Get the colormap from the provided name
        cmap = plt.cm.get_cmap(cmap_name, n)
        
        # Generate N colors as RGBA tuples
        colors = [cmap(i) for i in range(n)]
        
        return colors
    except ValueError:
        print(f"Error: '{cmap_name}' is not a valid colormap name.")
        return []
def copy_file_with_overwrite(src_file, dest_dir):
    # Ensure the destination directory exists
    os.makedirs(dest_dir, exist_ok=True)

    # Destination path for the file
    dest_path = os.path.join(dest_dir, os.path.basename(src_file))

    # Check if the destination file exists and is a symbolic link
    if os.path.islink(dest_path):
        print(f"Removing existing symbolic link: {dest_path}")
        os.unlink(dest_path)

    # Copy the file (overwrites if it exists)
    print(f"Copying {src_file} to {dest_path}")
    shutil.copy2(src_file, dest_path)
    # Verify copy and remove the source file
    if os.path.exists(dest_path):
        print(f"Successfully copied. Deleting source file: {src_file}")
        os.remove(src_file)
    else:
        print(f"Failed to copy {src_file} to {dest_path}. Source file not removed.")

def get_ew50_snr_fromcat(object_id=None, snr_threshold=5):
    """
    Function of get equivalent width (ew_50) and SNR of emission lines exceeding some threshold
    given grizli's initial concatenated output (Astropy Table), from all galaxies
    """
    idx2 = (np.where(grizli_cat['id'] == object_id))[0]
    lines = (grizli_cat['haslines'][idx2]).item()
    if isinstance(lines, bytes):  # If it's in bytes, decode
        lines = lines.decode('utf-8')
    
    # Now split into an array of line names
    lines_array = lines.split()
    
    # Extract SNR and EW50 values
    snr_ew50_list = []
    for line in lines_array:
        sn_col = f'sn_{line}'
        ew_col = f'ew50_{line}'
        
        if sn_col in grizli_cat.colnames and ew_col in grizli_cat.colnames:  # Ensure columns exist
            snr = float(grizli_cat[sn_col][idx2].item())  # Extract scalar value
            ew50 = float(grizli_cat[ew_col][idx2].item())  # Extract scalar value
            if snr >= snr_threshold:
                #print(line, snr, ew50)
                snr_ew50_list.append((line, snr, ew50))
    
    # Sort by SNR in descending order
    sorted_snr_ew50 = sorted(snr_ew50_list, key=lambda x: x[1], reverse=True)
    
    return sorted_snr_ew50
"""
Function to sort emission lines by SNR of lines
"""
def sorted_lines_by_snr(lines_prop_dicts=None):
    snr_dict = lines_prop_dicts['snr_line']

    # Sort emission lines by SNR (descending order)
    sorted_snr_list = sorted(snr_dict.items(), key=lambda x: x[1], reverse=True)
    
    # Construct the final sorted list of tuples with (Line_Name, SNR, EW50)
    sorted_snr_ew50 = [(line, snr, lines_prop_dicts['ew_50'][line]) for line, snr in sorted_snr_list]
    
    # Print the result
    return sorted_snr_ew50
def get_best_n_solutions_simple(zgrid=None, chi2=None, n=5, exclude_redshifts=[]):
    minima_indices, _ = scipy.signal.find_peaks(-1 * chi2, distance=50)
    minima_redshifts = zgrid[minima_indices]
    minima_chi_squared = chi2[minima_indices]

    # Exclude previously tried redshifts
    if exclude_redshifts:
        mask = ~np.isin(minima_redshifts, exclude_redshifts)
        minima_redshifts = minima_redshifts[mask]
        minima_chi_squared = minima_chi_squared[mask]

    # Sort minima by chi² value
    sorted_indices = np.argsort(minima_chi_squared)
    top_n_indices = sorted_indices[:n]

    # Get the n best solutions
    solution_redshifts = minima_redshifts[top_n_indices]
    solution_chi_squared = minima_chi_squared[top_n_indices]

    return solution_redshifts, solution_chi_squared



def iterative_redshift_search(object_id = 3763,zspec = 2.3994,showplot = False,do_find_bestredshift = True,
                              N_bestredshift = 4, max_iterations = 3,logprob_threshold=-6): 

    
    # Magnitude calculation
    idx = ((np.where(cat['id'] == object_id))[0])[0]
    
    try:
        mag_f200w = float(-2.5 * np.log10(cat['f200wn_tot_1'][idx]) + 23.9)
    except:
        #mag_f200w = float(-2.5 * np.log10(cat['f200wn_tot_2'][idx]) + 23.9)
        mag_f200w=np.nan
        pass

    try:
        mag_f150w = float(-2.5 * np.log10(cat['f150wn_tot_1'][idx]) + 23.9)
    except:
        mag_f150w=np.nan
        pass
    try:
        mag_f115w = float(-2.5 * np.log10(cat['f115wn_tot_1'][idx]) + 23.9)
    except:
        mag_f115w=np.nan
        pass
    
    if np.isfinite(mag_f200w):
        band = 'f200w'
        specific_mag = mag_f200w
    elif np.isfinite(mag_f150w):
        band = 'f150w'
        specific_mag = mag_f150w
    elif np.isfinite(mag_f115w):
        band='f115w'
        specific_mag=mag_f115w
    else:
        return -99 # none F115W, F150W and F200W are available
    # File paths
    full_file = f'{extract_dir}/{root}_{str(object_id).zfill(5)}.full.fits'
    mb_file = f'{extract_dir}/{root}_{str(object_id).zfill(5)}.beams.fits'
    args = np.load(fit_args_file, allow_pickle=True)[0]
    
    # Load FITS data
    if len(glob.glob(full_file)) == 0:
        return -99
    if len(glob.glob(full_file)) > 0:
        full_hdul = pyfits.open(full_file)
        zgrid = full_hdul['ZFIT_STACK'].data['zgrid']
        chi2 = full_hdul['ZFIT_STACK'].data['chi2'] / full_hdul['ZFIT_STACK'].header['DOF']
        
        # Main search loop
        attempt = 0
        excluded_redshifts = []
        z_best_new = -99
        
        while attempt < max_iterations and z_best_new == -99:
            print(f"\nIteration {attempt + 1}:")
            solution_redshifts, solution_chi_squared = get_best_n_solutions_simple(
                zgrid=zgrid, chi2=chi2, n=N_bestredshift, exclude_redshifts=excluded_redshifts
            )

            if (showplot) & (attempt==0):
                ## Plot chi2 distribution
                fig, ax = plt.subplots(figsize=(8, 6))
                
                # Plot the data
                ax.plot(zgrid, chi2, marker='none', linestyle='-', color='b', markersize=4, label=r'$\chi^2$ vs. Redshift',zorder=1)
                
                # Labels and title
                ax.set_xlabel('Redshift (zgrid)')
                ax.set_ylabel(r'$\chi^2$')
                ax.set_title(r'Redshift vs. $\chi^2$')
                
                # Grid and legend
                ax.grid(True, linestyle='--', alpha=0.6)
                
                #ax.set_xlim(np.min(solution_redshifts)-0.1,np.max(solution_redshifts)+0.1)
                if zspec is not None:
                    ax.axvline(zspec,color='crimson',ls=':')
                ax.set_ylim( np.min(chi2)-0.02, np.percentile(chi2,84))
                for k in range(len(solution_redshifts)):
                    ax.axvline(solution_redshifts[k],alpha=0.5,lw=2.5,zorder=-1,
                               color=line_colors[k],label='%.3f' % solution_redshifts[k])
        
                ax.legend(fontsize=myfontsize-5)
        
            print("Solution redshifts:", solution_redshifts)
            print("Solution chi-squared:", solution_chi_squared)
        
            logprob_arr = np.zeros(N_bestredshift) - 99
            z_map_arr = np.zeros(N_bestredshift) - 99
        
            for i, z_fit in enumerate(solution_redshifts):
                new_args = copy.deepcopy(args)
                new_args['zr'] = [z_fit]
                new_args['save_figures'] = False
                new_args['redshift_only'] = True
                new_args['write_fits_files'] = False
                new_args['verbose'] = False
        
                _mb, _, fit, tfit, line_hdu, cov_hdu = fitting.run_all(object_id, mb_files=[mb_file], **new_args)
                lines_prop_dicts = get_lines_fromcovar(cov_hdu=cov_hdu, min_snr=2, redshift=fit.meta["z_map"][0])
                sorted_snr_ew50 = sorted_lines_by_snr(lines_prop_dicts=lines_prop_dicts)
                ind=0 # the first line with highest SNR
                if sorted_snr_ew50:
                    logprob_arr[i] = compute_and_plot_joint_density(
                        fig=None, ax=None, target_line=sorted_snr_ew50[ind][0],
                        specific_ew=sorted_snr_ew50[ind][2], band=band,
                        specific_redshift=fit.meta["z_map"][0], specific_mag=specific_mag, showplot=False
                    )
                    z_map_arr[i] = fit.meta["z_map"][0]
                else:
                    print(f"No emission lines found for redshift {z_fit}.")
        
            # Check if any solution has a log likelihood > -8
            best_ind = np.where(logprob_arr > logprob_threshold)[0]
            if len(best_ind) > 0:
                z_map_arr_filtered = z_map_arr[best_ind]
                solution_chi_squared_filtered = solution_chi_squared[best_ind]
                z_best_new = z_map_arr_filtered[np.argmin(solution_chi_squared_filtered)]
                print(f"Valid solution found: z_best_new = {z_best_new}")
                break
        
            # Exclude the current batch of redshifts
            excluded_redshifts.extend(solution_redshifts)
            attempt += 1
        
        if z_best_new == -99:
            print("No valid redshift solution found after maximum iterations.")
        
        # Final fitting if a valid redshift is found
        if z_best_new > 0 and z_best_new != full_hdul['ZFIT_BEAM'].header['Z_MAP']:
            final_args = copy.deepcopy(args)
            final_args['zr'] = [z_best_new]
            final_args['verbose'] =False
            fitting.run_all(object_id, mb_files=[mb_file], **final_args)
        return z_best_new

def findobjects_2refit():
    object_id_toberefit = np.array([])
    for i in range(len(grizli_cat)):
        object_id=grizli_cat['id'][i]
        sorted_snr_ew50 = get_ew50_snr_fromcat(object_id=object_id, snr_threshold=5)
        if len(sorted_snr_ew50)>0 :
                   
            line_name =sorted_snr_ew50[0][0] # the first 0 mean the first line that has highest SNR,
            line_SNR = sorted_snr_ew50[0][1]
            line_EW =sorted_snr_ew50[0][2]
            if (line_name not in main_line_list) & ( line_EW>1):
                object_id_toberefit=np.hstack([object_id_toberefit ,object_id])
    np.savez('%s/%s_objectid_2refit.npz' % (new_redshift_dir,root),object_id=object_id_toberefit)
    
def main():
    npz=np.load('%s/%s_objectid_2refit.npz' % (new_redshift_dir,root))
    object_id_toberefit = npz['object_id']
    print(len(object_id_toberefit ))
    for object_id in object_id_toberefit:
        object_id=int(object_id)
        if object_id >0 :
            print(object_id)
            # check if the full.fits is still be symlink of initial grizli full.fits or the real file (not symlink)
            # if it is the real file then the redshift refinement has been done, do not need to do it again
            
            if (os.path.islink(new_redshift_dir+'/'+root+'_'+str(object_id).zfill(5) +'.full.fits')) | (len(glob.glob(new_redshift_dir+'/'+root+'_'+str(object_id).zfill(5) +'_nonewsolution.log'))==0):
                
                # extract the ew 50 of lines that have high enough SNR and EW (and sorted from high SNR to low)
                #sorted_snr_ew50 = get_ew50_snr_fromcat(object_id=object_id, snr_threshold=5)
                #print('object_id=',object_id)# ,sorted_snr_ew50)
                
                #if len(sorted_snr_ew50)>0 :
                   
                #line_name =sorted_snr_ew50[0][0] # the first 0 mean the first line that has highest SNR,
                #line_SNR = sorted_snr_ew50[0][1]
                #line_EW =sorted_snr_ew50[0][2]
                #if (line_name not in main_line_list) & ( line_EW>1):
                #try:   
                    #print(line_name,line_SNR,line_EW)
                try:
                    print('do  iterative_redshift_search')
                    z_best_new = iterative_redshift_search(object_id = object_id,zspec = None,
                                                           showplot = False,max_iterations = 2)
                except:
                    z_best_new=-99
                    pass
                #print('z_best_new',z_best_new)
                if z_best_new != -99:
                    new_file_list = glob.glob(f'*_{str(object_id).zfill(5)}*')
                    for file in new_file_list:
                        copy_file_with_overwrite(file, new_redshift_dir)
                if z_best_new ==-99:
                     with open(new_redshift_dir+'/'+root+'_'+str(object_id).zfill(5) +'_nonewsolution.log', "w") as file:
                         file.write("-99\n")



if __name__ == '__main__': 
    if len(glob.glob('%s/%s_objectid_2refit.npz' % (new_redshift_dir,root)))==0:
        findobjects_2refit()
    if len(glob.glob('%s/%s_objectid_2refit.npz' % (new_redshift_dir,root)))==1:
        main()             
    # Create fitting catalog
    consolidate_result=True
    
    if consolidate_result:
        print('combined all results')
        ids = Table.read(r'%s/%s-extracted.fits' % (new_redshift_dir,root))['NUMBER']
        results = vstack(
            [
                Table.read(r"%s/%s_%05d.row.fits" % (new_redshift_dir, root, i))
                for i in ids
            if os.path.exists(r"%s/%s_%05d.row.fits" % (new_redshift_dir,root, i))
            ]
            )

        results.rename_column('root', 'field')
        results.write(r"%s/%s_fitresults_useEWprior.fits" % (new_redshift_dir,root), overwrite=True)
