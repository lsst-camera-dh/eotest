from __future__ import print_function
from __future__ import absolute_import
import os
import sys
import subprocess
import numpy as np
import pylab
from .EOTestPlots import op_str


def latex_minus_value(value, error=None, format='%.2e'):
    """
    A latex entry for a single value, including optional error.
    """
    if value != value or error != error:
        # handle Nan's
        return "- \mbox{nan}"
    if value < 0:
        template = '+ \\num{' + format + '}'
    else:
        template = '- \\num{' + format + '}'
    result = template % np.abs(value)
    if error is not None:
        template = ' \pm \\num{' + format + '}'
        result += template % error
    return result


def _include_multipanel_png(pngfiles, frac_width=1.5, hspace=-1.9):
    return _include_png(pngfiles, frac_width=frac_width, hspace=hspace)


def _include_png(pngfiles, frac_width=1.3, hspace=-1):
    figure_template = """\\begin{figure}[H]
\\hspace{%.2fin}
%s
\\end{figure}
""" % (hspace, '%s')
    lines = []
    for item in pngfiles:
        if item.endswith('.png'):
            image_name = item[:-len('.png')]
        else:
            image_name = item
        lines.append('\\includegraphics[width=%s\\textwidth]{%s}'
                     % (frac_width, image_name))
    return figure_template % ('\\\\\n'.join(lines))


class EOTestReport(object):
    def __init__(self, eotest_plots, wl_dir, tex_file=None, qa_plot_files=None,
                 ccs_config_files=None, software_versions=None,
                 teststand_config=None, job_ids=None, sensor_grade_stats=None,
                 bnl_bias_offsets=None):
        self.plots = eotest_plots
        self.wl_dir = wl_dir
        if tex_file is None:
            # all files are written to the .tex without path, so write the .tex to the same place as images
            # then call pdflatex with the output_dir as the current working directory.
            self.tex_file = os.path.join(self.plots.output_dir, '%s_eotest_report.tex'%self.plots.sensor_id)
        else:
            self.tex_file = tex_file
        self.output = open(self.tex_file, 'w')
        self.qa_plot_files = qa_plot_files
        self.ccs_config_files = ccs_config_files
        self.software_versions = software_versions
        self.teststand_config = teststand_config
        self.job_ids = job_ids
        self.sensor_grade_stats = sensor_grade_stats
        self.bnl_bias_offsets = bnl_bias_offsets


    def make_figures(self):
        print("Creating eotest report figures...")
        funcs = ('fe55_dists',
                 'ptcs',
                 'gains',
                 'noise',
                 'full_well',
                 'linearity',
                 'linearity_resids',
                 'crosstalk_matrix',
                 'qe',
                 'psf_dists',
                 'persistence')
        for func in funcs:
            print("  %s" % func)
            try:
                exec('self.plots.%s()' % func)
                pylab.savefig('%s_%s.png' % (os.path.join(self.plots.output_dir, self.plots.sensor_id),
                                             func))
            except Exception as eobj:
                print("Error running %s():" % func)
                print("  ", eobj)
        self.plots.flat_fields(self.wl_dir)

    def _write_tex_preamble(self):
        self.output.write("""\documentclass{article}
\usepackage{graphicx}
\usepackage{amssymb}
\usepackage{siunitx}
\usepackage{float}
\usepackage[table,xcdraw]{xcolor}
\usepackage[margin=1in]{geometry}

\pagestyle{myheadings}

\\newcommand{\ok}{{\color[HTML]{009901} \checkmark}}
\\newcommand{\\fail}{{\color[HTML]{FE0000} $\\boldmath \\times$}}
\\newcommand{\electron}{{e$^-$}}
\\newcommand{\\twolinecell}[2][t]{%
  \\begin{tabular}[#1]{@{}l@{}}#2\end{tabular}}

\\begin{document}

""")

    def make_pdf(self):
        self._write_tex_preamble()
        sensor_id = self.plots.sensor_id
        #
        # Document title
        #
        self.output.write("\\title{Electro-Optical Test Results for Sensor %s}\n"
                          % sensor_id)
        self.output.write("\\maketitle\n")
        #
        # Write the summary table
        #
        self.output.write('\section{Summary}\n')
        self.output.write(self.plots.latex_table(hspace='-0.5in'))
        try:
            test_report_job_id = os.environ['LCATR_JOB_ID']
            self.output.write('Test Report Job ID: %s\n' % test_report_job_id)
        except KeyError:
            pass
        #
        # Write sensor grade stats
        #
        if self.sensor_grade_stats is not None:
            self.output.write('\n\section{Sensor Grade Statistics}\n')
            table = """\\begin{table}[h]
\\hspace{-0.5in}
\\begin{tabular}{""" + "|c"*len(self.sensor_grade_stats) + "|}\n"
            table += "\hline\n"
            table += ' & '.join(["\\texttt{%s}" % key for key in
                               self.sensor_grade_stats]) + "\\\\ \hline\n"
            fmts = ('%s', '%.2f', '%i', '%.2f', '%i', '%.2f', '%i',
                    '%.4f', '%i')
            table += ' & '.join([format % key for format, key in
                               zip(fmts, self.sensor_grade_stats.values())]) \
                               + "\\\\ \hline\n"
            table += '\end{tabular}\n\end{table}\n'
            self.output.write(table)
        #
        # Write BNL bias offsets
        #
        if self.bnl_bias_offsets is not None:
            # Format the table entries.
            entries = ['%i' % x for x in self.bnl_bias_offsets]

            # Compute the 8-channel medians and insert into table entries.
            entries.insert(0, '%.1f'  % np.median(self.bnl_bias_offsets[:8]))
            entries.insert(9, '%.1f'  % np.median(self.bnl_bias_offsets[8:]))

            # Write the table of bias value entries.
            table = """\\begin{table}[h]
\\centering
\\begin{tabular}{""" + "|c"*9 + "|}\n"
            table += "\hline\n"
            table += (' & '.join(('8-channel median',
                                "\\multicolumn{8}{c|}{BNL bias offset values, channels 1-8, 8-16}"))
                      + "\\\\ \hline\n")
            table += ' & '.join(entries[:9]) + "\\\\ \hline\n"
            table += ' & '.join(entries[9:]) + "\\\\ \hline\n"
            table += '\end{tabular}\n\end{table}\n'
            self.output.write(table)
        self.output.write('\\pagebreak\n\n')
        #
        # Read noise
        #
        self.output.write('\section{Read Noise}\n')
        self.output.write(self.plots.specs['CCD-007'].latex_table())
        self.output.write(_include_png(('%(sensor_id)s_noise' % locals(),)))
        self.output.write('\\pagebreak\n\n')
        #
        # Full Well and Nonlinearity
        #
        self.output.write('\section{Full Well and Nonlinearity}\n')
        self.output.write(self.plots.specs.latex_header() + '\n')
        self.output.write(self.plots.specs['CCD-008'].latex_entry() + '\n')
        self.output.write(self.plots.specs['CCD-009'].latex_entry() + '\n')
        self.output.write(self.plots.specs.latex_footer() + '\n')
        self.output.write("""\\begin{table}[!htbp]
\centering
\\begin{tabular}{|c|r|r|}
\hline
\\textbf{Amp} & \\textbf{Full Well} & \\twolinecell{\\textbf{Nonlinearity}\\\\(max. frac. dev.)} \\\\ \hline
""")
        full_well = self.plots.results['FULL_WELL']
        max_frac_dev = self.plots.results['MAX_FRAC_DEV']
        for amp in range(1, len(full_well)+1):
            my_full_well = op_str(full_well[amp-1], '%i')
            my_max_frac_dev = op_str(max_frac_dev[amp-1], '%.1e')
            self.output.write(" %i & $%s$ & $\\num{%s}$ \\\\ \hline\n"
                              % (amp, my_full_well, my_max_frac_dev))
        self.output.write("\\end{tabular}\n\\end{table}\n")

        """The following code is commented out because as of 2017/12/15 we are not producing
        the necessary data to produce these plots. It has not been *removed* as it is possible
        that this might be added in future. A config option to switch this is not appropraite
        as this is a 3rd party package, so this commenting out just lives on its own commit
        so that it can a) easily be undone, and b) won't be merged upsteam.

        self.output.write(_include_multipanel_png(('%(sensor_id)s_full_well'
                                                   % locals(),)))
        """

        self.output.write(_include_multipanel_png(('%(sensor_id)s_linearity'
                                                   % locals(),)))
        self.output.write(_include_multipanel_png(('%(sensor_id)s_linearity_resids'
                                                   % locals(),), frac_width=1.45,
                                                  hspace=-1.75))
        self.output.write('\\pagebreak\n\n')
        #
        # Mean Bias frame
        #
        mean_bias_file = '%(sensor_id)s_mean_bias' % locals()
        if os.path.isfile(mean_bias_file + '.png'):
            self.output.write('\section{Mean Bias Frame}\n')
            self.output.write(_include_png((mean_bias_file,)))
            self.output.write('\\pagebreak\n\n')
        #
        # Bright and Dark Defects mosaicked images.
        #
        bright_defects_file = '%(sensor_id)s_medianed_dark' % locals()
        dark_defects_file = '%(sensor_id)s_superflat_dark_defects' % locals()
        if (os.path.isfile(bright_defects_file + '.png') and
                os.path.isfile(dark_defects_file + '.png')):
            self.output.write('\section{Bright and Dark Defect Frames}')
            self.output.write(_include_png((bright_defects_file,)))
            self.output.write('\\pagebreak\n\n')
            self.output.write(_include_png((dark_defects_file,)))
            self.output.write('\\pagebreak\n\n')
        #
        # CTE
        #
        self.output.write('\section{Charge Transfer Efficiency}\n')
        self.output.write(self.plots.specs.latex_header() + '\n')
        self.output.write(self.plots.specs['CCD-010'].latex_entry() + '\n')
        self.output.write(self.plots.specs['CCD-011'].latex_entry() + '\n')
        self.output.write(self.plots.specs.latex_footer() + '\n')
        #
        # High flux level results
        #
        self.output.write('\subsection{High Flux}\n')
        self.output.write("""\\begin{table}[!htbp]
\centering
\\begin{tabular}{|c|l|l|}
\hline
\\textbf{Amp} & \\textbf{Serial CTE} & \\textbf{Parallel CTE} \\\\ \hline
""")
        scti = self.plots.results['CTI_HIGH_SERIAL']
        try:
            scti_err = self.plots.results['CTI_HIGH_SERIAL_ERROR']
        except KeyError:
            scti_err = np.zeros(len(scti))
        pcti = self.plots.results['CTI_HIGH_PARALLEL']
        try:
            pcti_err = self.plots.results['CTI_HIGH_PARALLEL_ERROR']
        except KeyError:
            pcti_err = np.zeros(len(scti))
        for amp in range(1, len(scti)+1):
            my_scti = latex_minus_value(scti[amp-1], error=scti_err[amp-1])
            my_pcti = latex_minus_value(pcti[amp-1], error=pcti_err[amp-1])
            self.output.write(" %(amp)i & $1%(my_scti)s$ & $1%(my_pcti)s$ \\\\ \hline\n" % locals())
        self.output.write("\\end{tabular}\n\\end{table}\n")
        #
        # Low flux level results
        #
        self.output.write('\subsection{Low Flux}\n')
        self.output.write("""\\begin{table}[!htbp]
\centering
\\begin{tabular}{|c|l|l|}
\hline
\\textbf{Amp} & \\textbf{Serial CTE} & \\textbf{Parallel CTE} \\\\ \hline
""")
        scti = self.plots.results['CTI_LOW_SERIAL']
        try:
            scti_err = self.plots.results['CTI_LOW_SERIAL_ERROR']
        except KeyError:
            scti_err = np.zeros(len(scti))
        pcti = self.plots.results['CTI_LOW_PARALLEL']
        try:
            pcti_err = self.plots.results['CTI_LOW_PARALLEL_ERROR']
        except KeyError:
            pcti_err = np.zeros(len(scti))
        for amp in range(1, len(scti)+1):
            my_scti = latex_minus_value(scti[amp-1], error=scti_err[amp-1])
            my_pcti = latex_minus_value(pcti[amp-1], error=pcti_err[amp-1])
            self.output.write(" %(amp)i & $1%(my_scti)s$ & $1%(my_pcti)s$ \\\\ \hline\n" % locals())
        self.output.write("\\end{tabular}\n\\end{table}\n")

        self.output.write('\\pagebreak\n\n')

        sflat_high = '%(sensor_id)s_superflat_high' % locals()
        sflat_low = '%(sensor_id)s_superflat_low' % locals()
        if (os.path.isfile(sflat_high + '.png') and
                os.path.isfile(sflat_low + '.png')):
            self.output.write(_include_png((sflat_high,)))
            self.output.write('\\pagebreak\n\n')
            self.output.write(_include_png((sflat_low,)))
            self.output.write('\\pagebreak\n\n')

        serial_profile_high = '%(sensor_id)s_serial_oscan_high' % locals()
        serial_profile_low = '%(sensor_id)s_serial_oscan_low' % locals()
        if (os.path.isfile(serial_profile_high + '.png') and
                os.path.isfile(serial_profile_low + '.png')):
            self.output.write(_include_png((serial_profile_high,)))
            self.output.write("""\\noindent
The blue points are the bias-subtracted column means; the red point is
the predicted signal in the columns used for the trailed charge based
on the measured CTI and the mean signal in the last column of the
imaging section; and the red dotted line is the signal based on the
CTE specification given the mean signal in the last imaging column.
""")
            self.output.write('\\pagebreak\n\n')
            self.output.write(_include_png((serial_profile_low,)))
            self.output.write("""\\noindent
The blue points are the bias-subtracted column means; the red point is
the predicted signal in the columns used for the trailed charge based
on the measured CTI and the mean signal in the last column of the
imaging section; and the red dotted line is the signal based on the
CTE specification given the mean signal in the last imaging column.
""")
            self.output.write('\\pagebreak\n\n')

        parallel_profile_high = '%(sensor_id)s_parallel_oscan_high' % locals()
        parallel_profile_low = '%(sensor_id)s_parallel_oscan_low' % locals()
        if (os.path.isfile(parallel_profile_high + '.png') and
                os.path.isfile(parallel_profile_low + '.png')):
            self.output.write(_include_png((parallel_profile_high,)))
            self.output.write("""\\noindent
The blue points are the bias-subtracted row means; the red point is
the predicted signal in the rows used for the trailed charge based on
the measured CTI and the mean signal in the last row of the imaging
section; and the red dotted line is the signal based on the CTE
specification given the mean signal in the last imaging row.
""")
            self.output.write('\\pagebreak\n\n')
            self.output.write(_include_png((parallel_profile_low,)))
            self.output.write("""\\noindent
The blue points are the bias-subtracted row means; the red point is
the predicted signal in the rows used for the trailed charge based on
the measured CTI and the mean signal in the last row of the imaging
section; and the red dotted line is the signal based on the CTE
specification given the mean signal in the last imaging row.
""")
            self.output.write('\\pagebreak\n\n')
        #
        # Crosstalk
        #
        self.output.write('\section{Crosstalk}\n')
        self.output.write(self.plots.specs.latex_header() + '\n')
        self.output.write(self.plots.specs['CCD-013'].latex_entry() + '\n')
        self.output.write(self.plots.specs.latex_footer() + '\n')
        crosstalk_image = '%(sensor_id)s_crosstalk_matrix' % locals()
        if os.path.isfile(crosstalk_image + '.png'):
            self.output.write(_include_png((crosstalk_image,)))
        self.output.write('\\pagebreak\n\n')

        """The following code is commented out because as of 2017/12/15 we are not producing
        the necessary data to produce these plots. It has not been *removed* as it is possible
        that this might be added in future. A config option to switch this is not appropraite
        as this is a 3rd party package, so this commenting out just lives on its own commit
        so that it can a) easily be undone, and b) won't be merged upsteam.
        #
        # QE
        #
        self.output.write('\section{Quantum Efficiency}\n')
        self.output.write(self.plots.specs.latex_header() + '\n')
        self.output.write(self.plots.specs['CCD-021'].latex_entry() + '\n')
        self.output.write(self.plots.specs['CCD-022'].latex_entry() + '\n')
        self.output.write(self.plots.specs['CCD-023'].latex_entry() + '\n')
        self.output.write(self.plots.specs['CCD-024'].latex_entry() + '\n')
        self.output.write(self.plots.specs['CCD-025'].latex_entry() + '\n')
        self.output.write(self.plots.specs['CCD-026'].latex_entry() + '\n')
        self.output.write(self.plots.specs.latex_footer() + '\n')
        self.output.write(_include_png(('%(sensor_id)s_qe' % locals(),)))
        self.output.write('\\pagebreak\n\n')
        """
        #
        # PRNU
        #
        self.output.write('\section{Photoresponse Non-uniformity}\n')
        self.output.write(self.plots.specs['CCD-027'].latex_table())
        flats = []
        for wl in self.plots.prnu_wls:
            flat = '%(sensor_id)s_%(wl)04inm_flat' % locals()
            if os.path.isfile(flat + '.png'):
                self.output.write(self.plots.specs.prnu_specs[wl].latex_table())
                self.output.write(_include_png((flat,)))
                self.output.write('\\pagebreak\n\n')
        #
        # Point Spread Function
        #
        self.output.write('\section{Point Spread Function}\n')
        self.output.write(self.plots.specs['CCD-028'].latex_table())
        self.output.write(_include_multipanel_png(('%(sensor_id)s_psf_dists'
                                                   % locals(),)))
        self.output.write('\\pagebreak\n\n')
        #
        # Image persistence plots
        #
        persistence_image = '%(sensor_id)s_persistence' % locals()
        if os.path.isfile(persistence_image + '.png'):
            self.output.write('\section{Image Persistence}\n')
            self.output.write(_include_multipanel_png((persistence_image,)))
            self.output.write('\\pagebreak\n\n')
        #
        # Fe55 gains and PTC
        #
        self.output.write('\section{System Gain and Photon Transfer Curves}\n')
        self.output.write(_include_multipanel_png(('%(sensor_id)s_fe55_dists'
                                                   % locals(),)))
        self.output.write(_include_png(('%(sensor_id)s_gains' % locals(),)))
        self.output.write('\\pagebreak\n\n')
        if os.path.isfile('%(sensor_id)s_ptcs.png' % locals()):
            self.output.write(_include_multipanel_png(('%(sensor_id)s_ptcs'
                                                       % locals(),)))
        self.output.write('\\pagebreak\n\n')
        #
        # Fe55 zoom
        #
        fe55_zoom = '%(sensor_id)s_fe55_zoom' % locals()
        if os.path.isfile(fe55_zoom + '.png'):
            self.output.write('\section{Fe55 zoom on segment 1}\n')
            self.output.write(_include_png((fe55_zoom,)))
            self.output.write('\\pagebreak\n\n')
        #
        # Fe55 aperture flux vs pixel number.
        #
        fe55_apflux_serial = '%(sensor_id)s_fe55_apflux_serial' % locals()
        fe55_apflux_parallel = '%(sensor_id)s_fe55_apflux_parallel' % locals()
        if (os.path.isfile(fe55_apflux_serial + '.png') and
                os.path.isfile(fe55_apflux_parallel + '.png')):
            self.output.write('\section{Fe55 aperture flux vs pixel number}\n')
            self.output.write(_include_png((fe55_apflux_serial,)))
            self.output.write(_include_png((fe55_apflux_parallel,)))
            self.output.write('\\pagebreak\n\n')
        #
        # Fe55 p3-p5 statistics.
        #
        fe55_p3_p5_hists = '%(sensor_id)s_fe55_p3_p5_hists' % locals()
        fe55_p3_p5_profiles = '%(sensor_id)s_fe55_p3_p5_profiles' % locals()
        if (os.path.isfile(fe55_p3_p5_hists + '.png') and
                os.path.isfile(fe55_p3_p5_profiles + '.png')):
            self.output.write('\section{Fe55 p3-p5 statistics}\n')
            self.output.write(_include_png((fe55_p3_p5_hists,)))
            self.output.write(_include_png((fe55_p3_p5_profiles,)))
            self.output.write('\\pagebreak\n\n')
        #
        # QA plots
        #
        if self.qa_plot_files is not None:
            self.output.write('\section{QA plots}\n')
            for item in self.qa_plot_files:
                self.output.write(_include_png((item,), frac_width=1.1,
                                               hspace=-0.5))
                self.output.write('\\pagebreak\n\n')
        #
        # eTraveler ActivityIDs
        if self.job_ids is not None:
            self.output.write('\section{eTraveler Activity IDs}\n')
            self.output.write("""\\begin{table}[!htbp]
\centering
\\begin{tabular}{|l|r|}
\hline
\\textbf{Job Name} & \\textbf{activityId} \\\\ \hline
""")
            values = np.array([int(x) for x in list(self.job_ids.values())])
            keys = list(self.job_ids.keys())
            index = np.argsort(values)
            print("sorted job id indexes:", index)
            for i in index:
                job_name = keys[i].replace('_', '\_')
                self.output.write("%s & %i \\\\ \hline\n"
                                  % (job_name, values[i]))
            self.output.write("\\end{tabular}\n\\end{table}\n")
            self.output.write('\\pagebreak\n\n')
        #
        # Software versions
        #
        if self.software_versions is not None:
            self.output.write('\section{Software Versions}\n')
            self.output.write('\\begin{description}\n')
            for key, value in list(self.software_versions.items()):
                self.output.write('\\item[%s] %s\n' % (key.replace('_', '\_'),
                                                       value.replace('_', '\_')))
            self.output.write('\end{description}\n')
        #
        # Test stand configuration
        #
        if self.teststand_config is not None:
            self.output.write('\section{Test Stand Configuration}\n')
            self.output.write('\\begin{description}\n')
            for key, value in list(self.teststand_config.items()):
                self.output.write('\\item[%s] %s\n' %
                                  (key.replace('_', '\_'),
                                   str(value).replace('_', '\_')))
            self.output.write('\end{description}\n')
#        #
#        # CCS configuration files
#        #
#        if self.ccs_config_files is not None:
#            self.output.write('\section{CCS Configurations}\n')
#            for key, full_path in self.ccs_config_files.items():
#                sys.stdout.flush()
#                keyname = key.replace('_', '\_')
#                filename = os.path.basename(full_path).replace('_', '\_')
#                self.output.write('\subsection{%s: %s}\n' % (keyname, filename))
#                self.output.write('\\begin{verbatim}\n')
#                for line in open(full_path):
#                    self.output.write(line)
#                self.output.write('\end{verbatim}\n')
#                self.output.write('\\pagebreak\n\n')
        #
        # End and close the document.
        #
        self.output.write("\\end{document}\n")
        self.output.close()

        subprocess.call('pdflatex %s' % self.tex_file, shell=True, cwd=os.path.dirname(self.tex_file))


if __name__ == '__main__':
    from .EOTestPlots import EOTestPlots
    plots = EOTestPlots('000-00', 'results', output_dir='plots')
    report = EOTestReport(plots, './sensorData/000-00/lambda/debug')
    report.make_figures()
    report.make_pdf()
