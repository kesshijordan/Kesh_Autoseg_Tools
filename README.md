# Automatic Segmentation Pipeline for Virtual Dissection of White Matter Fascicle Models

This pipeline shows an example implementation of a tractography segmentation pipeline that combines two complementary approaches: shape recognition and end point constraint with a region-of-interest (ROI) approach.

Streamlines included in a fascicle model are identified based on
1. anatomical connectivity priors based on regions of interest (ROIs) in the subject’s native space
2. shape priors based on 3D streamline bundle atlases applied using [Recobundles](https://pubmed.ncbi.nlm.nih.gov/28712994/) implemented in [DIPY](https://dipy.org)

Note that this approach can be tuned for individual bundles by modifying the inclusion/exclusion ROIs used, the bundle templates, and the Recobundles parameters.

A related methodology (Jordan, 2018 unpublished) to create a custom bundle template has also been provided to complement this pipeline ([Python Notebook Link](make_custom_streamline_template.ipynb)).

## Contents of Repository

### ["Run Pipeline" Python File](command_line_apply_pipeline_example.py)

This is the main file that runs the pipeline. The code is written in python, using the open source [Diffusion Imaging in Python](https://dipy.org/) package and assumes that
1.  all pre-processing up to the Freesurfer Parcellation has already been run, aligned in the subject’s diffusion space and stored as nifti-format ROIs and
2. A whole-brain tractography streamline dataset has been run and is stored in .trk format in the subject’s native diffusion space.

This code depends on the following publicly available python packages: os, nibabel, numpy, dipy, and pandas. The bundle template choice and the rules for ROI targeting are specified through two csv files so that no python skills are needed to create a new bundle segmentation or modify an existing one after the pipeline is set up to run.

The command-line interface requires 3 arguments: the whole brain tractography streamline dataset of the case to be segmented (case_whole_brain), the freesurfer parcellation aparc aseg file in nift format (case_aparc_aseg), and the path to the folder where the output will be stored (case_putpath).

### [Bundle Atlas CSV File](bundle_template_lookup.csv)
This file ([bundle_template_lookup.csv](bundle_template_lookup.csv)) is used to define the bundle atlases that will be used to constrain the bundle extraction via [Recobundles](https://pubmed.ncbi.nlm.nih.gov/28712994/) shape recognition. The pipeline will attempt to segment any tract that has a bundle template in this file, so it is also the way to specify which bundles to segment when running this pipeline. The path to the bundle template needs to be specified, which also serves as a selection of that bundle for segmentation.

#### Publicly Available Bundle Templates
- [Subset of (Yeh, 2018) HCP Tractography Atlas bundles in MNI Space](https://figshare.com/articles/Simple_model_bundle_atlas_for_RecoBundles/6483614)
- [Full set of (Yeh, 2018) HCP Tractography Atlas bundles](http://brain.labsolver.org/diffusion-mri-templates/tractography)

### [ROI Target Rules CSV File](ROI_target_rules_lookup.csv)
This file ([ROI_target_rules_lookup.csv](ROI_target_rules_lookup.csv)) is used to define the targeting rules that constrain the endpoints of the bundle segmentation. The positive numbers indicate inclusion ROIs with numbers indicating sets of Freesurfer ROIs that are combined with an “OR” operation. Distinct enumerated sets of ROIs are used for targeting with an “AND” operation . The negative value indicates an exclusion ROI, which are used for targeting with a “NOT” operation. For each column, streamlines must intersect each positively enumerated group of ROIs, and not intersect with any -1 ROIs.  


### [Custom Streamline Template Demo](make_custom_streamline_template.ipynb)
This notebook ([make_custom_streamline_template.ipynb](make_custom_streamline_template.ipynb)) demonstrates how a streamline clustering method [Quickbundles](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3518823/), a bundle-based registration method [Streamline Linear Registration](https://www.sciencedirect.com/science/article/abs/pii/S1053811915003961), and an outlier rejection method [CCI](https://pubmed.ncbi.nlm.nih.gov/28940825/) can be used to define a custom bundle template from a set of fascicle models (e.g. manually curated). This is an unpublished methodology presented at a hackathon (Jordan, 2018).

## Adapting this code for your own use
Make sure that you update all of the paths in the python code & csv's to reflect your own setup. The load/save and streamline manipulation functions in this example are specific to manually curated bundles that were stored prior to recent standardization of formats. These functions may need to be replaced to reflect the specifics of your bundles.

## Adding or Removing a Bundle
To remove a bundle from the segmentation, delete the path from the file_name column in the bundle atlas csv. To add a new bundle, add a row to the bundle template csv with a bundle name and path and add a column to the freesurfer csv file specifying the ROI targeting rules using exactly the same name for the new bundle in both csv files.

## Custom Bundle Templates
If you are interested in generating your own bundle templates, code is available in the [Custom Bundle Template python notebook](make_custom_streamline_template.ipynb). This validation used publicly available bundle templates but the creation of custom templates is an avenue to explore in tailoring the methodology to the population under investigation.

## Application to Other Cohorts
The specific implementation provided can be easily adapted for a variety of use cases by modifying the ROI priors using the ROI Target Rules csv and the shape priors with a different bundle template. The creation of custom templates using the [make_custom_streamline_template.ipynb](make_custom_streamline_template.ipynb) notebook can also be used to tailor a segmentation to a cohort that does not fit the adult healthy control models well.

## References

Garyfallidis, Eleftherios (2018): Simple model bundle atlas for RecoBundles. figshare. Dataset. https://doi.org/10.6084/m9.figshare.6483614.v1

Yeh, F. C., S. Panesar, D. Fernandes, A. Meola, M. Yoshino, J. C. Fernandez-Miranda, J. M. Vettel, and T. Verstynen. "Population-averaged atlas of the macroscale human structural connectome and its network topology." NeuroImage (2018).

Jordan KM. RecoBundles Templates. Presented at BrainHack Global at Indiana University; May 4, 2018; Bloomington.

Garyfallidis E, Ocegueda O, Wassermann D, et al. Robust and efficient linear registration of white-matter fascicles in the space of streamlines. NeuroImage 2015;117:124–40.

Garyfallidis E, Côté M-A, Rheault F, et al. Recognition of white matter bundles using local and global streamline-based registration and clustering. NeuroImage 2018;170:283–95.

Garyfallidis E, Brett M, Correia MM, Williams GB, Nimmo-Smith I. QuickBundles, a method for tractography simplification. Front Neurosci 2012;6:175.
