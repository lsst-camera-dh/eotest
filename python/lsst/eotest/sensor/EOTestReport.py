import os
import sys
import subprocess
import numpy as np
import pylab

def latex_minus_value(value, error=None, format='%.2e'):
    """
    A latex entry for a single value, including optional error.
    """
    if value < 0:
        template = '+ \\num{' + format + '}'
    else:
        template = '- \\num{' + format + '}'
    result = template % np.abs(value)
    if error is not None:
        template = ' \pm \\num{' + format + '}'
        result += template % error
    return result

def _include_png(pngfiles, frac_height=0.6):
    figure_template = """\\begin{figure}[H]
\\centering
%s
\\end{figure}
"""
    lines = []
    for item in pngfiles:
        if item.endswith('.png'):
            image_name = item[:-len('.png')]
        else:
            image_name = item
        lines.append('\\includegraphics[height=%s\\textheight]{%s}' 
                     % (frac_height, image_name))
    return figure_template % ('\\\\\n'.join(lines))

class EOTestReport(object):
    def __init__(self, eotest_plots, wl_dir, tex_file=None, qa_plot_files=None,
                 ccs_config_files=None, software_versions=None,
                 teststand_config=None):
        self.plots = eotest_plots
        self.wl_dir = wl_dir
        if tex_file is None:
            self.tex_file = '%s_eotest_report.tex' % self.plots.sensor_id
        else:
            self.tex_file = tex_file
        self.output = open(self.tex_file, 'w')
        self.qa_plot_files = qa_plot_files
        self.ccs_config_files = ccs_config_files
        self.software_versions = software_versions
        self.teststand_config = teststand_config
    def make_figures(self):
        print "Creating eotest report figures..."
        funcs = ('fe55_dists',
                 'ptcs',
                 'gains',
                 'noise',
                 'full_well',
                 'linearity',
                 'crosstalk_matrix',
                 'qe',
                 'psf_dists',
                 'persistence')
        for func in funcs:
            print "  %s" % func
            try:
                exec('self.plots.%s()' % func)
                pylab.savefig('%s_%s.png' % (self.plots.sensor_id, func))
            except Exception, eobj:
                print "Error running %s():" % func
                print "  ", eobj
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
        self.output.write(self.plots.latex_table())
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
        self.output.write(self.plots.specs.latex_footer()+ '\n')
        self.output.write("""\\begin{table}[!htbp]
\centering
\\begin{tabular}{|c|r|r|}
\hline
\\textbf{Amp} & \\textbf{Full Well} & \\twolinecell{\\textbf{Nonlinearity}\\\\(max. frac. dev.)} \\\\ \hline
""")
        full_well = self.plots.results['FULL_WELL']
        max_frac_dev = self.plots.results['MAX_FRAC_DEV']
        for amp in range(1, 17):
            my_full_well = full_well[amp-1]
            my_max_frac_dev = max_frac_dev[amp-1]
            self.output.write(" %i & $%i$ & $\\num{%.1e}$ \\\\ \hline\n" 
                              % (amp, my_full_well, my_max_frac_dev))
        self.output.write("\\end{tabular}\n\\end{table}\n")
        self.output.write(_include_png(('%(sensor_id)s_full_well' % locals(),)))
        self.output.write(_include_png(('%(sensor_id)s_linearity' % locals(),)))
        self.output.write('\\pagebreak\n\n')
        #
        # CTE
        #
        self.output.write('\section{Charge Transfer Efficiency}\n')
        self.output.write(self.plots.specs.latex_header() + '\n')
        self.output.write(self.plots.specs['CCD-010'].latex_entry() + '\n')
        self.output.write(self.plots.specs['CCD-011'].latex_entry() + '\n')
        self.output.write(self.plots.specs.latex_footer()+ '\n')
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
            scti_err = np.zeros(16)
        pcti = self.plots.results['CTI_HIGH_PARALLEL']
        try:
            pcti_err = self.plots.results['CTI_HIGH_PARALLEL_ERROR']
        except KeyError:
            pcti_err = np.zeros(16)
        for amp in range(1, 17):
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
            scti_err = np.zeros(16)
        pcti = self.plots.results['CTI_LOW_PARALLEL']
        try:
            pcti_err = self.plots.results['CTI_LOW_PARALLEL_ERROR']
        except KeyError:
            pcti_err = np.zeros(16)
        for amp in range(1, 17):
            my_scti = latex_minus_value(scti[amp-1], error=scti_err[amp-1])
            my_pcti = latex_minus_value(pcti[amp-1], error=pcti_err[amp-1])
            self.output.write(" %(amp)i & $1%(my_scti)s$ & $1%(my_pcti)s$ \\\\ \hline\n" % locals())
        self.output.write("\\end{tabular}\n\\end{table}\n")

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
        self.output.write(_include_png(('%(sensor_id)s_psf_dists' % locals(),)))
        self.output.write('\\pagebreak\n\n')
        #
        # Image persistence plots
        #
        persistence_image = '%(sensor_id)s_persistence' % locals()
        if os.path.isfile(persistence_image + '.png'):
            self.output.write('\section{Image Persistence}\n')
            self.output.write(_include_png((persistence_image,)))
            self.output.write('\\pagebreak\n\n')
        #
        # Fe55 gains and PTC
        #
        self.output.write('\section{System Gain and Photon Transfer Curves}\n')
        self.output.write(_include_png(('%(sensor_id)s_fe55_dists' 
                                        % locals(),)))
        self.output.write(_include_png(('%(sensor_id)s_gains' % locals(),),
                                       frac_height=0.45))
        self.output.write('\\pagebreak\n\n')
        if os.path.isfile('%(sensor_id)s_ptcs.png' % locals()):
            self.output.write(_include_png(('%(sensor_id)s_ptcs' % locals(),)))
        self.output.write('\\pagebreak\n\n')
        #
        # QA plots
        #
        if self.qa_plot_files is not None:
            self.output.write('\section{QA plots}\n')
            for item in self.qa_plot_files:
                self.output.write(_include_png((item,), frac_height=0.9))
                self.output.write('\\pagebreak\n\n')
        #
        # Software versions
        #
        if self.software_versions is not None:
            self.output.write('\section{Software Versions}\n')
            self.output.write('\\begin{description}\n')
            for key, value in self.software_versions.items():
                self.output.write('\\item[%s] %s\n' % (key.replace('_', '\_'),
                                                       value.replace('_', '\_')))
            self.output.write('\end{description}\n')
        #
        # Test stand configuration
        #
        if self.teststand_config is not None:
            self.output.write('\section{Test Stand Configuration}\n')
            self.output.write('\\begin{description}\n')
            for key, value in self.teststand_config.items():
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

        subprocess.call('pdflatex %s' % self.tex_file, shell=True)

if __name__ == '__main__':
    from EOTestPlots import EOTestPlots
    plots = EOTestPlots('000-00', 'results', output_dir='plots')
    report = EOTestReport(plots, './sensorData/000-00/lambda/debug')
    report.make_figures()
    report.make_pdf()
