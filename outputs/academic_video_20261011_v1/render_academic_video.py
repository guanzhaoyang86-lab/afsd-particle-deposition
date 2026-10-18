"""Minimal academic 3D movie of a fixed, verified chain of saved AFSD states."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import subprocess
import time

os.environ.setdefault('NUMBA_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from render_kernels import base_scene, spheres

HERE = Path(__file__).resolve().parent
FFMPEG = Path('C:/Program Files/JianyingPro/Apps/5.7.0.11527/ffmpeg.exe')
W, H, FPS, SECONDS = 1920, 1080, 30, 26
INK, MUTED, LINE = (37, 50, 63), (106, 119, 131), (223, 229, 234)
BLUE, GREY, TAN, ORANGE = (66, 128, 163), (190, 201, 208), (178, 165, 140), (217, 112, 47)
FONT, BOLD = 'C:/Windows/Fonts/msyh.ttc', 'C:/Windows/Fonts/msyhbd.ttc'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def smooth(v):
    v = float(np.clip(v, 0, 1))
    return v*v*(3-2*v)


def basis(azimuth, elevation):
    a, e = math.radians(azimuth), math.radians(elevation)
    view = np.array([math.cos(e)*math.cos(a), math.cos(e)*math.sin(a), math.sin(e)])
    right = np.cross([0., 0., 1.], view)
    right /= np.linalg.norm(right)
    return right, np.cross(view, right), view


class Movie:
    def __init__(self):
        with np.load(HERE/'data_frames.npz', allow_pickle=False) as z:
            self.d = {k: z[k].copy() for k in z.files}
        self.m = json.loads((HERE/'data_manifest.json').read_text(encoding='utf-8'))
        self.xyz = self.d['xyz_m']*1000
        self.origins = self.d['tool_origin_m']*1000
        self.ids = self.d['particle_id']
        self.times = self.d['time_s']
        self.feed = self.d['feed_mask'].astype(bool)
        self.eligible = self.d['eligible_mask'].astype(bool)
        self.crossed = self.d['supplied_crossed_mask'].astype(bool)
        self.under = self.d['under_shoulder_mask'].astype(bool)
        self.params = self.m['parameters']
        self.radius = self.params['tool_radius']*1000
        self.bore = self.params['bore_width']*500
        self.tool_height = (self.params['tool_top']-self.params['layer'])*1000
        self.fonts = {}
        self.frame_map = []
        assert self.xyz.shape == (len(self.times), len(self.ids), 3)
        assert np.all(np.diff(self.times) > 0) and np.isfinite(self.xyz).all()
        self.input_hashes = {p.name: sha(p) for p in [HERE/'data_frames.npz', HERE/'data_manifest.json']}

    def font(self, size, bold=False):
        key = size, bold
        if key not in self.fonts:
            self.fonts[key] = ImageFont.truetype(BOLD if bold else FONT, size)
        return self.fonts[key]

    def text(self, draw, xy, text, size=23, color=INK, bold=False, anchor=None):
        draw.text(xy, text, font=self.font(size, bold), fill=color, anchor=anchor)

    def source_index(self, n):
        u = float(np.clip((n/FPS-4)/16, 0, 1))
        target = self.times[0]+u*(self.times[-1]-self.times[0])
        return min(len(self.times)-1, max(0, int(np.searchsorted(self.times, target+1e-14, side='right')-1)))

    def colors(self, k):
        c = np.empty((len(self.ids), 3), dtype=float)
        c[:] = GREY
        c[self.feed & ~self.eligible] = TAN
        c[self.eligible] = BLUE
        c[self.crossed[k]] = ORANGE
        c[self.under[k]] = (57, 143, 116)
        return c

    @staticmethod
    def project(point, target, right, up, scale, w, h):
        p = np.asarray(point)-target
        return np.array([w/2+np.dot(p, right)*scale, h/2-np.dot(p, up)*scale])

    def view(self, k, n, detail=False):
        sec = n/FPS
        az = -69+13*smooth(sec/4)+8*smooth((sec-20)/6)
        elev = 21 if not detail else 18
        right, up, view = basis(az, elev)
        width, height = (1050, 710) if not detail else (704, 584)
        target = np.array([1.5, 0., 2.15]) if not detail else self.origins[k]+np.array([0., 0., .04])
        scale = 63. if not detail else 196.
        origin, rotation = self.origins[k], self.d['tool_rotation'][k]
        if detail:
            rgb = np.full((height, width, 3), 255, dtype=np.uint8)
            depth = np.full((height, width), 1e9, dtype=np.float32)
            q = self.xyz[k]-origin
            mask = (q[:, 2] >= -.66) & (q[:, 2] <= .72) & (np.linalg.norm(q[:, :2], axis=1) < 1.95)
        else:
            # Display-only half section of the rigid tool. Material remains intact.
            rgb, depth = base_scene(width, height, scale, target, right, up, view,
                                    origin, rotation, self.radius, self.bore, self.tool_height, 0.)
            mask = np.ones(len(self.ids), dtype=bool)
        spheres(rgb, depth, np.ascontiguousarray(self.xyz[k, mask]), np.ascontiguousarray(self.colors(k)[mask]),
                scale, target, right, up, view, .105 if detail else .185)
        im = Image.fromarray(rgb)
        draw = ImageDraw.Draw(im)
        proj = lambda p: self.project(p, target, right, up, scale, width, height)
        if detail:
            b = self.bore
            square = np.array([[-b,-b,0.],[b,-b,0.],[b,b,0.],[-b,b,0.],[-b,-b,0.]])@rotation.T+origin
            points = [tuple(proj(p)) for p in square]
            # Dashed reference plane outline is an annotation, not a material surface.
            for a, b in zip(points[:-1], points[1:]):
                a, b = np.array(a), np.array(b)
                length = np.linalg.norm(b-a)
                for s in np.arange(0., length, 13.):
                    start = a+(b-a)*(s/length)
                    end = a+(b-a)*(min(s+7, length)/length)
                    draw.line((tuple(start),tuple(end)), fill=(135,153,166), width=2)
            # All four historical crossing centres remain visible as annotated rings.
            marked=np.flatnonzero(self.crossed[k])
            median_y=float(np.median([proj(self.xyz[k,i])[1] for i in marked])) if len(marked) else height/2
            for i in marked:
                x, y = proj(self.xyz[k,i])
                draw.ellipse((x-11,y-11,x+11,y+11), outline=ORANGE, width=3)
                sx = -1 if x < width/2 else 1
                dx, dy = sx*56, 48 if y >= median_y else -48
                draw.line((x+sx*12,y, x+dx,y+dy), fill=ORANGE, width=1)
                label=str(int(self.ids[i])); lx,ly=x+dx+sx*5,y+dy-10
                anchor='ra' if sx<0 else None
                box=draw.textbbox((lx,ly),label,font=self.font(18),anchor=anchor)
                draw.rectangle((box[0]-3,box[1]-2,box[2]+3,box[3]+2),fill='white')
                self.text(draw,(lx,ly),label,18,ORANGE,anchor=anchor)
            self.text(draw, (18, 9), '工具隐藏 · 孔口下方局部', 20, MUTED)
            self.text(draw, (18, 40), '虚线：方孔出口平面', 19, MUTED)
            self.text(draw, (width-18,height-34), '圆环：已穿孔粒子中心', 19, MUTED, anchor='ra')
        else:
            def callout(text, xy, point):
                self.text(draw, xy, text, 21, MUTED)
                tw=draw.textlength(text,font=self.font(21))
                p=proj(point)
                a=(xy[0]+tw+12,xy[1]+16)
                draw.line((a,(a[0]+24,a[1]),tuple(p)), fill=(160,174,185), width=1)
                draw.ellipse((p[0]-2,p[1]-2,p[0]+2,p[1]+2),fill=MUTED)
            callout('送料棒', (52,35), [-.1,0.,5.65])
            callout('旋转工具 / 剖视', (52,236), [-1.7,.25,2.35])
            callout('可变形基板', (52,430), [-2.75,-2.5,-.7])
            # Orthographic coordinate triad with world directions.
            x0,y0=100,645
            for axis,label,color in zip(np.eye(3),['x','y','z'],[MUTED]*3):
                dx,dy=np.dot(axis,right)*42,-np.dot(axis,up)*42
                draw.line((x0,y0,x0+dx,y0+dy), fill=color,width=2)
                self.text(draw,(x0+dx*1.25,y0+dy*1.25),label,19,color,anchor='mm')
        bar_mm = 1. if not detail else .5
        xx,yy=(890,673) if not detail else (52,548)
        ll=bar_mm*scale
        draw.line((xx,yy,xx+ll,yy),fill=INK,width=3)
        for x in (xx,xx+ll): draw.line((x,yy-5,x,yy+5),fill=INK,width=2)
        self.text(draw,(xx+ll/2,yy+10),f'{bar_mm:g} mm',17,MUTED,anchor='ma')
        return im

    def frame(self, n, record=True):
        k=self.source_index(n)
        im=Image.new('RGB',(W,H),'white')
        draw=ImageDraw.Draw(im)
        self.text(draw,(66,37),'AFSD  /  3D PARTICLE OBSERVATION',21,MUTED)
        self.text(draw,(64,72),'材料流动与孔口穿越',40,INK,True)
        self.text(draw,(1854,38),'模拟时间 / SIMULATION TIME',19,MUTED,anchor='ra')
        self.text(draw,(1854,70),f'{self.times[k]*1000:.4f} ms',39,INK,True,anchor='ra')
        draw.line((65,139,1855,139),fill=LINE,width=2)
        self.text(draw,(66,159),'a   整体三维结构',24,INK,True)
        self.text(draw,(1160,159),'b   孔口局部放大',24,INK,True)
        im.paste(self.view(k,n),(52,196))
        im.paste(self.view(k,n,True),(1152,214))
        draw=ImageDraw.Draw(im)
        draw.line((1120,186,1120,894),fill=LINE,width=1)
        self.text(draw,(1168,823),'累计新料穿孔',22,MUTED)
        self.text(draw,(1360,810),str(int(self.crossed[k].sum())),48,ORANGE,True)
        self.text(draw,(1450,823),'肩下占据',22,MUTED)
        self.text(draw,(1630,810),str(int(self.under[k].sum())),48,INK,True)
        self.text(draw,(1168,879),'本时间窗内：穿孔数量未增加，尚未形成肩下沉积。',20,MUTED)
        labels=[('基板',GREY),('初始下端棒料',TAN),('待穿孔新料',BLUE),('已穿孔新料',ORANGE)]
        x=78
        for label,c in labels:
            draw.ellipse((x,921,x+14,935),fill=c)
            self.text(draw,(x+24,913),label,20,MUTED)
            x+=draw.textlength(label,font=self.font(20))+74
        self.text(draw,(1855,915),f'{len(self.ids):,} 粒子  /  {len(self.times)} 个真实保存状态',20,MUTED,anchor='ra')
        left,right,yy=67,1855,981
        draw.line((left,yy,right,yy),fill=LINE,width=3)
        u=(self.times[k]-self.times[0])/(self.times[-1]-self.times[0])
        if u>0: draw.line((left,yy,left+(right-left)*u,yy),fill=BLUE,width=3)
        phase='初始保存状态 · 相机绕视' if n/FPS<4 else '真实状态慢放' if n/FPS<20 else '末态观察 · 相机绕视'
        self.text(draw,(67,949),phase,19,MUTED)
        self.text(draw,(1855,949),f'{self.times[0]*1000:.3f} — {self.times[-1]*1000:.3f} ms',19,MUTED,anchor='ra')
        self.text(draw,(67,1012),'工程近似  ·  真实保存位置  ·  位移未放大  ·  球体为粒子标记',20,MUTED)
        self.text(draw,(1855,1012),'右图圆环为透视定位标注  /  原比例',19,MUTED,anchor='ra')
        if record:
            self.frame_map.append(dict(video_frame=n,video_time_s=n/FPS,source_frame=k,
                physical_time_s=float(self.times[k]),crossed_count=int(self.crossed[k].sum()),
                under_shoulder_count=int(self.under[k].sum())))
        return im


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--preview',action='store_true'); args=ap.parse_args()
    movie=Movie(); began=time.perf_counter()
    if args.preview:
        samples=[0,360,779]
        sheet=Image.new('RGB',(1440,810),'white')
        for j,n in enumerate(samples):
            frame=movie.frame(n,False); frame.save(HERE/f'preview_{n:04d}.png')
            sheet.paste(frame.resize((720,405),Image.Resampling.LANCZOS),((j%2)*720,(j//2)*405))
        sheet.save(HERE/'storyboard.jpg',quality=94)
        print('PREVIEWS_READY',flush=True); return
    video=HERE/'AFSD_学术三维观察_20261011.mp4'
    command=[str(FFMPEG),'-y','-hide_banner','-loglevel','warning','-f','rawvideo','-pix_fmt','rgb24',
        '-s',f'{W}x{H}','-r',str(FPS),'-i','-','-an','-c:v','h264_nvenc','-preset','p4','-rc','vbr','-cq','19','-b:v','8M',
        '-threads','2','-pix_fmt','yuv420p','-movflags','+faststart',str(video)]
    with (HERE/'encoder.log').open('wb') as log:
        proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=log,stderr=log,
            creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        try:
            for n in range(FPS*SECONDS):
                frame=movie.frame(n); proc.stdin.write(frame.tobytes())
                if n in (0,360,779): frame.save(HERE/f'frame_{n:04d}.png')
                if n%90==0: print(f'RENDER {n}/{FPS*SECONDS} elapsed={time.perf_counter()-began:.1f}s',flush=True)
            proc.stdin.close(); code=proc.wait(timeout=60)
            assert code==0,(code,(HERE/'encoder.log').read_text(errors='replace'))
        except BaseException:
            proc.kill(); proc.wait(); raise
    assert all(sha(HERE/p)==h for p,h in movie.input_hashes.items())
    (HERE/'frame_sources.json').write_text(json.dumps(movie.frame_map,indent=2),encoding='utf-8')
    receipt=dict(status='RENDERED_REAL_SAVED_STATES_PENDING_DECODE_QA',physical_solver_evaluations=0,
        input_sha256=movie.input_hashes,renderer_sha256=sha(Path(__file__)),kernel_sha256=sha(HERE/'render_kernels.py'),
        source_count=len(movie.times),source_particle_count=len(movie.ids),source_time_s=[float(movie.times[0]),float(movie.times[-1])],
        physical_position_interpolation=False,displacement_magnification=1,
        glyph_radius_mm=dict(overall=.185,detail=.105),detail_glyphs_reduced_for_visibility=True,
        tool_display='world-y half section in overview; hidden in local inset; geometry and poses from saved states',
        particle_visibility='overview all particles; local saved centres in declared spatial crop; crossing rings are x-ray annotations',
        local_crop='relative outlet z in [-0.66,0.72] mm; cylindrical radius < 1.95 mm',
        playback='4 s first state camera orbit; 16 s uniform requested physical time with previous saved sample held; 6 s final state camera orbit',
        video=dict(path=str(video),sha256=sha(video),bytes=video.stat().st_size,width=W,height=H,fps=FPS,frames=FPS*SECONDS,duration_s=SECONDS),
        encode_command=command,elapsed_wall_seconds=time.perf_counter()-began)
    (HERE/'render_receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(receipt['video'],ensure_ascii=False),flush=True)


if __name__=='__main__': main()
