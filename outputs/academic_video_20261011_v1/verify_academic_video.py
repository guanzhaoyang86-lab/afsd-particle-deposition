"""Decode the exported movie and check its saved-state frame mapping."""
from pathlib import Path
import json
import subprocess
import numpy as np
from PIL import Image
from render_academic_video import HERE, FFMPEG, Movie, sha

receipt_path=HERE/'render_receipt.json'
receipt=json.loads(receipt_path.read_text(encoding='utf-8'))
video=Path(receipt['video']['path'])
qa=HERE/'qa'; qa.mkdir(exist_ok=True)
flags=getattr(subprocess,'CREATE_NO_WINDOW',0)
r=subprocess.run([str(FFMPEG),'-hide_banner','-nostats','-i',str(video),'-map','0:v:0',
    '-progress','pipe:1','-f','null','NUL'],capture_output=True,text=True,creationflags=flags,timeout=60)
(qa/'decode.log').write_text(r.stderr+'\n'+r.stdout,encoding='utf-8')
assert r.returncode==0,r.stderr
for expected in ['1920x1080','30 fps','Duration: 00:00:26.00','Video: h264','yuv420p']:
    assert expected in r.stderr,expected
assert 'frame=780' in r.stdout and 'progress=end' in r.stdout
indices=[0,360,779]
select='select='+'+'.join(f'eq(n\\,{n})' for n in indices)
r=subprocess.run([str(FFMPEG),'-y','-v','error','-i',str(video),'-vf',select,'-vsync','0',
    str(qa/'decoded_%02d.png')],capture_output=True,text=True,creationflags=flags,timeout=60)
assert r.returncode==0,r.stderr
movie=Movie()
mapping=json.loads((HERE/'frame_sources.json').read_text(encoding='utf-8'))
assert len(mapping)==780
prior=-1
for n,row in enumerate(mapping):
    k=row['source_frame']
    assert row['video_frame']==n and k>=prior and k==movie.source_index(n)
    assert row['physical_time_s']==movie.times[k]
    assert row['crossed_count']==int(movie.crossed[k].sum())
    assert row['under_shoulder_count']==int(movie.under[k].sum())
    prior=k
assert mapping[0]['source_frame']==0 and mapping[-1]['source_frame']==len(movie.times)-1
assert sha(HERE/'data_frames.npz')==movie.m['npz_sha256']
assert sha(video)==receipt['video']['sha256']
assert all(sha(HERE/p)==h for p,h in receipt['input_sha256'].items())
assert sha(HERE/'render_academic_video.py')==receipt['renderer_sha256']
assert sha(HERE/'render_kernels.py')==receipt['kernel_sha256']
errors=[]
contact=Image.new('RGB',(1440,810),'white')
for j,n in enumerate(indices):
    decoded=Image.open(qa/f'decoded_{j+1:02d}.png').convert('RGB')
    source=Image.open(HERE/f'frame_{n:04d}.png').convert('RGB')
    error=float(np.abs(np.asarray(decoded).astype(float)-np.asarray(source).astype(float)).mean())
    assert error<4,(n,error)
    errors.append(dict(frame=n,mean_absolute_rgb_error=error))
    contact.paste(decoded.resize((720,405),Image.Resampling.LANCZOS),((j%2)*720,(j//2)*405))
contact.save(qa/'decoded_contact_sheet.jpg',quality=95)
result=dict(status='PASS_DECODE_AND_SAVED_STATE_MAPPING',decode_exit_code=0,frame_count=780,duration_s=26,
    width=1920,height=1080,fps=30,codec='H264',pixel_format='yuv420p',
    source_frame_count=len(movie.times),displayed_distinct_source_frames=len(set(x['source_frame'] for x in mapping)),
    crossed_counts=sorted(set(x['crossed_count'] for x in mapping)),
    under_shoulder_counts=sorted(set(x['under_shoulder_count'] for x in mapping)),
    all_frame_mappings_checked=True,physical_position_interpolation=False,sources_unchanged=True,
    decoded_sample_errors=errors,visual_review='PENDING',video_sha256=sha(video))
(qa/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
