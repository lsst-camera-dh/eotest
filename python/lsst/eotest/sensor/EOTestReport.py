import os
import subprocess
import pylab
from EOTestPlots import EOTestPlots

def _include_png(pngfiles, frac_height=0.6):
    figure_template = """\\begin{figure}[!htbp]
\\centering
%s
\\end{figure}
"""
    lines = []
    for item in pngfiles:
        lines.append('\\includegraphics[height=%s\\textheight]{%s}' 
                     % (frac_height, item.strip('.png')))
    return figure_template % ('\\\\\n'.join(lines))

class EOTestReport(object):
    def __init__(self, eotest_plots, wl_dir, tex_file=None):
        self.plots = eotest_plots
        self.wl_dir = wl_dir
        if tex_file is None:
            self.tex_file = '%s_eotest_report.tex' % self.plots.sensor_id
        else:
            self.tex_file = tex_file
        self.output = open(tex_file, 'w')
    def make_figures(self):
        print "Creating eotest report figures..."
        funcs = ('fe55_dists',
                 'ptcs',
                 'gains',
                 'noise',
                 'linearity',
                 'crosstalk_matrix',
                 'qe',
                 'psf_dists')
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
\usepackage[table,xcdraw]{xcolor}
\usepackage[margin=1in]{geometry}

\pagestyle{myheadings}

\\begin{document}

\\newcommand{\ok}{{\color[HTML]{009901} \checkmark}}
\\newcommand{\\fail}{{\color[HTML]{FE0000} $\\boldmath \\times$}}
\\newcommand{\electron}{{e$^-$}}
\\newcommand{\\twolinecell}[2][t]{%
  \\begin{tabular}[#1]{@{}l@{}}#2\end{tabular}}
""")
    def make_pdf(self):
        self._write_tex_preamble()
        sensor_id = self.plots.sensor_id
        #
        # Write the summary table
        #
        self.output.write('\section{Summary}\n')
        self.output.write(plots.latex_table())
        self.output.write('\\pagebreak\n')
        #
        # Read noise
        #
        self.output.write('\section{Read Noise}\n')
        self.output.write(plots.specs['CCD-007'].latex_table())
        self.output.write(_include_png(('%(sensor_id)s_noise' % locals(),)))
        self.output.write('\\pagebreak\n')
        #
        # Full Well and Nonlinearity
        #
        self.output.write('\section{Full Well and Nonlinearity}\n')
        self.output.write(plots.specs.latex_header() + '\n')
        self.output.write(plots.specs['CCD-008'].latex_entry() + '\n')
        self.output.write(plots.specs['CCD-009'].latex_entry() + '\n')
        self.output.write(plots.specs.latex_footer()+ '\n')
        self.output.write(_include_png(('%(sensor_id)s_linearity' % locals(),)))
        self.output.write('\\pagebreak\n')
        #
        # Crosstalk
        #
        self.output.write('\section{Crosstalk}\n')
        self.output.write(plots.specs.latex_header() + '\n')
        self.output.write(plots.specs['CCD-013'].latex_entry() + '\n')
        self.output.write(plots.specs.latex_footer() + '\n')
        self.output.write(_include_png(('%(sensor_id)s_crosstalk_matrix' 
                                  % locals(),)))
        self.output.write('\\pagebreak\n')
        #
        # QE
        #
        self.output.write('\section{Quantum Efficiency}\n')
        self.output.write(plots.specs.latex_header() + '\n')
        self.output.write(plots.specs['CCD-021'].latex_entry() + '\n')
        self.output.write(plots.specs['CCD-022'].latex_entry() + '\n')
        self.output.write(plots.specs['CCD-023'].latex_entry() + '\n')
        self.output.write(plots.specs['CCD-024'].latex_entry() + '\n')
        self.output.write(plots.specs['CCD-025'].latex_entry() + '\n')
        self.output.write(plots.specs['CCD-026'].latex_entry() + '\n')
        self.output.write(plots.specs.latex_footer() + '\n')
        self.output.write(_include_png(('%(sensor_id)s_qe' % locals(),)))
        self.output.write('\\pagebreak\n')
        #
        # PRNU
        #
        self.output.write('\section{Photresponse Non-uniformity}\n')
        self.output.write(plots.specs['CCD-027'].latex_table())
        flats = []
        for wl in plots.prnu_wls:
            flat = '%(sensor_id)s_%(wl)04inm_flat' % locals()
            if os.path.isfile(flat + '.png'):
                self.output.write(_include_png((flat,)))
        self.output.write('\\pagebreak\n')
        #
        # Point Spread Function
        #
        self.output.write('\section{Point Spread Function}\n')
        self.output.write(plots.specs['CCD-028'].latex_table())
        self.output.write(_include_png(('%(sensor_id)s_psf_dists' % locals(),)))
        self.output.write('\\pagebreak\n')
        #
        # Fe55 gains and PTC
        #
        self.output.write('\section{System Gain and Photon Transfer Curves}\n')
        self.output.write(_include_png(('%(sensor_id)s_fe55_dists' 
                                        % locals(),)))
        self.output.write(_include_png(('%(sensor_id)s_gains' % locals(),),
                                       frac_height=0.45))
        self.output.write('\\pagebreak\n')
        self.output.write(_include_png(('%(sensor_id)s_ptcs' % locals(),)))
        
        self.output.write("\\end{document}\n")
        self.output.close()

        subprocess.call('pdflatex %s' % self.tex_file, shell=True)

if __name__ == '__main__':
    plots = EOTestPlots('000-00', 'results', output_dir='plots')
    report = EOTestReport(plots, './sensorData/000-00/lambda/debug')
    report.make_figures()
    report.make_pdf()
