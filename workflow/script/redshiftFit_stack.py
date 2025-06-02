#! /usr/bin/env python

# Import packages
import os
import warnings
import argparse
from multiprocessing import Pool
from astropy.table import Table, vstack
import numpy as np
# Import grizli
import grizli
from grizli import fitting
from grizli.pipeline import auto_script
import grizli.pipeline.photoz
import pickle
# Silence warnings
warnings.filterwarnings('ignore')
import scipy
print(scipy.__version__)

def main():
    # Parse arguements
    parser = argparse.ArgumentParser()
    parser.add_argument('fieldname', type=str)
    parser.add_argument('--ncpu', type=int, default=1)
    parser.add_argument('--part',type=int)
    args = parser.parse_args()
    fname = args.fieldname
    ncpu = args.ncpu
    chunk_index = args.part
    # Print version and step
    print(f'Fitting {fname}')
    print(f'grizli:{grizli.__version__}')

    # Get paths and get fields
    main = os.getcwd()
    fields = os.path.join(main, 'FIELDS')
    home = os.path.join(fields, fname)

    # Subdirectories
    # plots = os.path.join(home, 'Plots')
    extract = os.path.join(home, 'Extractions_stack')
    os.chdir(extract)

    # Return if directory is empty
    if not os.path.exists(f'{fname}-extracted.fits'):
        print('No extracted spectra found')
        return

    # Generate fit arguements
    pline = {
        'kernel': 'square',
        'pixfrac': 0.5,
        'pixscale': 0.04,
        'size': 8,
        'wcs': None,
    }
    #args=auto_script.generate_fit_params(
    #        fit_only_beams=True,
    #        fit_beams=True,
    #        run_fit=True,
    #        poly_order=7,
    #        sys_err=0.1,
    #        fcontam=1,
    #        zr=[0.1,12],
    #        fit_trace_shift=True,
    #        save_file='fit_args.npy',
    #    pline=pline, field_root=fname, min_sens=0.0, min_mask=0.0
    #)
    args = np.load('fit_args.npy', allow_pickle=True)[0]
    import glob
    if len(glob.glob('outthere-niriss-%s_eazyphot.pkl' % (fname)))==0:
        aper_ix = 1
        total_flux_column = 'flux'
        get_external_photometry = False
        int_ez = grizli.pipeline.photoz.eazy_photoz(fname, force=False, object_only=True,
                      apply_background=True, aper_ix=aper_ix,
                      apply_prior=True, beta_prior=True,
                      get_external_photometry=get_external_photometry,
                      external_limits=3, external_sys_err=0.3, external_timeout=300,
                      sys_err=args['sys_err'], z_step=0.01, z_min=0.01, z_max=12,
                      total_flux=total_flux_column)
        t0 = args['t0']
        int_phot_obj = grizli.pipeline.photoz.EazyPhot(int_ez, grizli_templates=t0, zgrid=int_ez.zgrid)
        with open('outthere-niriss-%s_eazyphot.pkl' % (fname), "wb") as file:  # Open the file in binary write mode
            pickle.dump(int_phot_obj, file)
    if len(glob.glob('outthere-niriss-%s_eazyphot.pkl' % (fname)))==1:
        with open('outthere-niriss-%s_eazyphot.pkl' % (fname),'rb') as f:
            int_phot_obj = pickle.load(f)
    order=0
    #args_withPhot=auto_script.generate_fit_params(
    #    fit_only_beams=True,
    #    fit_beams=True,
    #    run_fit=True,
    #    poly_order=7,
    #    sys_err=0.1,
    #    fcontam=1,
    #    zr=[0.1,12],
    #    fit_trace_shift=True,
    #    use_phot_obj=True,
    #    scale_photometry=order+1,
    #    phot_obj=int_phot_obj,
    #    save_file='fit_args_wNIRISSphot.npy',
    #    pline=pline, field_root=fname, min_sens=0.00, min_mask=0.0
    #)
    args_withPhot = np.load('fit_args_wNIRISSphot.npy', allow_pickle=True)[0]
    #args_withPhot['only_stacks']=True
    #args_withPhot['fit_stacks']=True
    #args_withPhot['fit_beams']=False
    #np.save('fit_args_wNIRISSphot_fitstack.npy', [args_withPhot])
    # Get IDs
    _ids = Table.read(f'{fname}-extracted.fits')['NUMBER']
    # Split into 15 chunks
    chunks = np.array_split(_ids, 15)
    ids = chunks[chunk_index] # from 0 to 14
    files = glob.glob(extract+'/*full*.fits')
    done_id = [int(os.path.basename(files[i]).replace(fname+'_','').replace('.full.fits','')) for i in range(len(files))]
    notexist_id = list(set(ids)-set(done_id))
    notexist_id = np.sort(np.array(notexist_id))
    print('notexist_id',notexist_id)
    if len(notexist_id)>0:
        # Multiprocessing pool
        with Pool(processes=ncpu) as pool:
        #    for result in pool.imap_unordered(zfit,notexist_id, chunksize=4):
        #        pass
            
            #result = pool.map_async(zfit, ids)
            result = pool.map_async(zfit, notexist_id)
            _ = result.get()
            pool.close()
            pool.join()
    #for object_id in notexist_id:
    #    zfit(object_id)
    # Create fitting catalog
    #results = vstack(
    #    [
    #        Table.read(f'{fname}_{str(i).zfill(5)}.row.fits')
    #        for i in ids
    #        if os.path.exists(f'{fname}_{str(i).zfill(5)}.row.fits')
    #    ]
    #)
    #results.rename_column('root', 'field')
    #results.write(f'{fname}_fitresults.fits', overwrite=True)


# Fit Redshift
def zfit(i):
    # Fit
    fitting.run_all_parallel(i, args_file='fit_args_wNIRISSphot.npy', verbose=True,fit_only_beams=False, only_stacks= False,fit_beams=True)
           

    # # Create oned spectrum figure
    # fig = mb.oned_figure()
    # fig.savefig(os.path.join(extract_plots,f'{i}_oned.png'))
    # pyplot.close(fig)

    # # Create 2d spectrum figure
    # hdu,fig = mb.drizzle_grisms_and_PAs(size=38, scale=0.5, diff=True, kernel='square', pixfrac=0.5)
    # fig.savefig(os.path.join(extract_plots,f'{i}_twod.png'))
    # pyplot.close(fig)


if __name__ == '__main__':
    #print("Hello")
    main()
