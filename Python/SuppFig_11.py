#!/usr/bin/env python3
# Documentation (Supp Fig 11):
# This script runs all 14 samples through DAPI segmentation, protein background correction, cell
# classification and cross-K analysis of CD8A+/Granzyme B+ T cells and CD274+ macrophages.
# The helper modules are in Functions/ (figs11_segmentation_helpers.py, figs11_analysis_helpers.py)
# and the final threshold JSON is AuxFiles/figs11_thresholds.json. It produces per-sample outputs and
# combined result tables. Sample information and reported statistics are provided in Table S20.
"""Run Figure S11 for every sample in AuxFiles/figs11_thresholds.json.

Run: python SuppFig_11.py --run-id FigS11_run1

segmentation -> corrected intensities -> cell gates -> cross-K.
The two helper modules in the Functions/ directory contain the processing functions.
"""

from pathlib import Path
import argparse
import sys

code_dir = "/projects/Pregnancy_placenta/Code_Publication/Python/"  # root of the code hierarchy
HERE = Path(code_dir)
sys.path.insert(0, str(HERE / 'Functions'))


# SETTINGS: shared by every sample. The JSON supplies all sample-specific gates.
# -----------------------------------------------------------------------------
THRESHOLDS_FILE = HERE / 'AuxFiles' / 'figs11_thresholds.json'
CELLPOSE_VERSION = '4.0.6'
ARCSINH_COFACTOR = 5.0
ZSCORE_DDOF = 1
SEED = 20260930
N_PERMUTATIONS = 9999
RADII_UM = list(range(10, 151, 5))

# Usually empty. Only needed to choose among acquisitions of the same sample.
# Values are paths relative to the Pregnancy_complications_Miltenyi project.
IMAGE_OVERRIDES = {}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id', default='FigS11_run1', help='New output run name.')
    parser.add_argument('--segmentation-run', default=None,
                        help='Optional existing segmentation run to reuse for every sample.')
    args = parser.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import figs11_segmentation_helpers as segmentation
    import figs11_analysis_helpers as analysis

    plt.rcParams.update({'font.family': 'Arial', 'pdf.fonttype': 42, 'ps.fonttype': 42})
    segmentation.require_compute()

    # 1. LOAD ALL SAMPLES: no TMA/ROI/ROW editing or separate sample submissions.
    samples = analysis.load_cohort(THRESHOLDS_FILE)
    sample_ids = [sample['sample_id'] for sample in samples]
    output_folder = segmentation.result_directory(args.run_id, sample_ids[0]).parent
    segmentation_run = args.segmentation_run or f'{args.run_id}_segmentation'
    if output_folder.exists():
        raise FileExistsError(f'Use a new run ID; results already exist: {output_folder}')
    if args.segmentation_run is None:
        segmentation_folder = segmentation.result_directory(segmentation_run, sample_ids[0]).parent
        if segmentation_folder.exists():
            raise FileExistsError(f'Segmentation output already exists: {segmentation_folder}')
    core_tables, radius_tables = [], []
    print(f'Processing {len(samples)} samples. Results: {output_folder}', flush=True)

    # Process one sample at a time to keep memory use bounded; all are automatic.
    for position, sample in enumerate(samples, start=1):
        sample_id = sample['sample_id']
        tma, roi, row = sample['TMA'], sample['ROI'], sample['ROW']
        image_override = IMAGE_OVERRIDES.get(sample_id)
        print(f'\n[{position}/{len(samples)}] {sample_id}', flush=True)

        # 2. SEGMENT DAPI AND EXPAND NUCLEI BY 5 µm.
        if args.segmentation_run is None:
            print('  Segmenting the registered DAPI image...', flush=True)
            image = segmentation.selected_image(tma, roi, image_relative=image_override)
            pixel_size = segmentation.get_pixel_size_um(image)
            segmentation.run_segmentation(
                tma, roi, row, run_id=segmentation_run,
                pixel_size_um=pixel_size, expected_cellpose_version=CELLPOSE_VERSION,
                image_relative=image_override,
            )

        # 3. APPLY SAVED PIXEL THRESHOLDS AND MEASURE CORRECTED CELL INTENSITIES.
        print('  Measuring corrected protein intensities...', flush=True)
        inputs = analysis.load_sample(
            tma, roi, row, segmentation_run=segmentation_run,
            thresholds_path=THRESHOLDS_FILE, image_relative=image_override,
        )
        cells = analysis.measure_corrected_cells(inputs)
        del inputs

        # 4. APPLY SAVED Z-SCORE GATES AND ASSIGN CELL POPULATIONS.
        print('  Applying cell gates...', flush=True)
        cells = analysis.classify_cells(cells, cofactor=ARCSINH_COFACTOR, ddof=ZSCORE_DDOF)

        # 5. RUN THE FIGURE S11 CROSS-K COMPARISON.
        print('  Running cross-K permutations...', flush=True)
        results = analysis.analyze_crossk(
            cells, seed=SEED, n_permutations=N_PERMUTATIONS, radii_um=RADII_UM,
        )

        # 6. SAVE THIS SAMPLE'S FINAL CELLS, FIGURES AND COMPLETE SOURCE DATA.
        figures = {
            'marker_gates': analysis.plot_gate_distributions(cells),
            'cell_populations': analysis.plot_cell_populations(cells),
            'crossK': analysis.plot_crossk(results),
        }
        saved_folder = analysis.save_results(cells, results, run_id=args.run_id, figures=figures)
        core_tables.append(results['core_stats'])
        radius_tables.append(results['radius_stats'])
        for figure in figures.values():
            plt.close(figure)
        del cells, results, figures
        print(f'  Saved: {saved_folder}', flush=True)

    # 7. SAVE COMBINED TABLES FOR ALL 14 SAMPLES, IN THE JSON'S ORIGINAL ORDER.
    summary = analysis.save_cohort_results(
        core_tables, radius_tables, output_folder=output_folder, sample_ids=sample_ids,
    )
    segmentation.write_json(output_folder / 'cohort_run.json', dict(
        samples=sample_ids, segmentation_run=segmentation_run,
        seed=SEED, n_permutations=N_PERMUTATIONS, radii_um=RADII_UM,
        code=segmentation.helper_fingerprints(), execution=segmentation.execution_record(),
    ))
    print('\nCompleted all samples.', flush=True)
    print(summary[['sample_id', 'p_global', 'max_T']].to_string(index=False))
    print(f'Combined tables: {output_folder}', flush=True)
