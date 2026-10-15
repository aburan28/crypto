"""Render native-derived per-session values; perform no statistical analysis."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

root=Path(__file__).resolve().parents[1]
data=json.loads((root/'results/session-summaries.json').read_text())
pdfmetrics.registerFont(TTFont('DejaVu','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
styles=getSampleStyleSheet()
for style in styles.byName.values(): style.fontName='DejaVu'
styles['BodyText'].fontSize=9;styles['BodyText'].leading=13
for session in data['sessions']:
    directory=root/'reports'/session['session'];directory.mkdir(parents=True,exist_ok=True)
    cells=session['cells']
    fig,ax=plt.subplots(figsize=(7,3.5))
    ax.bar([c['arm'] for c in cells],[c['mean_recorded_process_cpu_seconds'] for c in cells],color='#258491')
    ax.set_ylabel('Mean recorded process user + system CPU seconds')
    ax.set_title(session['session']+' · all recorded attempt statuses')
    fig.text(.05,.015,'L0; descriptive CPU. Native work incompletely priced; no IC/rho ratio.',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,1))
    fig.savefig(directory/'process-cpu.svg',bbox_inches='tight')
    fig.savefig(directory/'process-cpu.png',dpi=150,bbox_inches='tight');plt.close(fig)
    text=[f"# {session['session']}",'',session['scope'],'',
          'L0, co-scheduled. No complete common-price IC/rho claim.',
          '',f'[Native-derived summary](../../../../research/conductor_gap_cm_orbits_20261015/results/session-summaries.json) · [Raw session](../../../../research/conductor_gap_cm_orbits_20261015/sessions/{session["session"]}/records.jsonl)',
          '',f"Records SHA-256: `{session['records_sha256']}`",f"Executable SHA-256: `{session['binary_sha256']}`",'']
    story=[Paragraph(session['session'],styles['Title']),Paragraph(session['scope'],styles['BodyText']),Spacer(1,10),
           Image(str(directory/'process-cpu.png'),width=480,height=240)]
    for cell in cells:
        line=f"{cell['arm']}: {cell['attempts']} attempts; statuses {json.dumps(cell['statuses'],sort_keys=True)}. Mean recorded process CPU {cell['mean_recorded_process_cpu_seconds']:.6f} s."
        text.extend([line,'']);story.extend([Paragraph(line,styles['BodyText']),Spacer(1,8)])
        for phase in cell['online_phases']:
            line=f"{phase['phase']}: {phase['mean_wall_seconds']:.6f} s mean wall, {phase['observations']} available online receipts."
            text.append(line);story.append(Paragraph(line,styles['BodyText']))
    text.extend(['','The graph contains all attempt statuses; phase means omit missing receipts and name their observation counts. Failures are retained. No successful-only speed ratio is computed.'])
    (directory/'REPORT.md').write_text('\n'.join(text)+'\n')
    story.extend([Spacer(1,12),Paragraph('L0, co-scheduled; incomplete native price. Failures and censored attempts are retained. These CPU and phase diagnostics do not establish ECDLP hardness or an admitted IC/rho speedup.',styles['BodyText'])])
    SimpleDocTemplate(str(directory/'REPORT.pdf')).build(story)
