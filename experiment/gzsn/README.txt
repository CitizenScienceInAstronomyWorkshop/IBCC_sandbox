Currently-working scripts:

runGZSN_strongprior - Use this one for doing current experiments. Use the settings file to adjust subsampling.
runGZSN_pkgmain2

Obsolete-looking files:

runGZSN
runGZSN_flatprior
runGZSN_nosubsampling
runGZSN_oneFold
runGZSN_pkgmain


Things to do:
0. Ditch the subsampling if IBCCVB is still reasonable. Or show a simple subsampling scheme as a comparison. Then present IBCCDiff as a way around this.
1. Can we interpret the ground truth so that the mean is no more than 70% reliable?
    a) Look at labelling and subsampling currently used.
    b) Run with no subsampling or minimal subsampling.
2. Does the number of -ve examples skew the ROC? 
3. Look at the subsampling script. Can we simplify it to get an IBCC improvement over mean? E.g. only use workers with n classifications. 
4. Run GZSN with simplified subsampling.
5. Run with no subsampling against IBCCDiff. Should get improvement without any subsampling.
 