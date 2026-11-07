"""Akdeniz Tool for HU Calculation, version 1.0.0.

Interactive measurement of Hounsfield units (HU) in a 3D region of interest (ROI)
placed on a spinal CT series (DICOM folder).

Workflow: open a DICOM folder, place eight ROI points on sagittal slices, press
Enter to compute the ROI statistics (copied to the clipboard), or press P to
export the summary figure for publication.

Developed by Muhammad Ateya, Department of Neurosurgery, Akdeniz University.
"""
import sys
import os
import numpy as np
import pydicom
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from PIL import Image, ImageTk
import matplotlib
matplotlib.use('TkAgg')  # Tk backend must be selected explicitly (needed for PyInstaller builds)
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.path import Path
from scipy.stats import norm
from mpl_toolkits.mplot3d import Axes3D  # registers the 3D projection in older matplotlib versions

# Disable matplotlib's default toolbar shortcuts, which conflict with the tool's
# own keys (C, V, Z, N, P, U, Backspace). Keymaps absent from the installed
# matplotlib version are skipped.
for _keymap in ('keymap.back', 'keymap.forward', 'keymap.zoom', 'keymap.pan',
                'keymap.save', 'keymap.grid', 'keymap.grid_minor',
                'keymap.yscale', 'keymap.xscale', 'keymap.all_axes'):
    if _keymap in plt.rcParams:
        plt.rcParams[_keymap] = []

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# =============================================================================
# Publication figure (key P, after all eight ROI points are placed)
#   Panel A: sagittal and coronal CT slices with the 3D ROI outline.
#   Panel B: HU distribution of the ROI voxels with mean, median and FWHM mean
#            HU; the Gaussian-based FWHM interval is shaded.
# Export settings: white background, Arial (Liberation Sans or DejaVu Sans if
# Arial is unavailable), text >= 8 pt, 300 dpi. No DICOM header fields are read
# or drawn, so no patient identifier appears in the figure.
# Colors are from the Okabe-Ito colorblind-safe palette; line styles also differ.
# =============================================================================
PUB_DPI = 300
PUB_WINDOW = (-200, 1000)          # CT display window (HU), panel A
PUB_OUTLINE = '#E69F00'            # 3D ROI outline
PUB_MEAN, PUB_MEDIAN, PUB_FWHM_MEAN = '#D55E00', '#0072B2', '#009E73'
PUB_SHADE, PUB_BARS = '#56B4E9', '#999999'


def make_publication_figure(volume, mask, voxel_spacing, data, out_path, dpi=PUB_DPI):
    """Save the two-panel ROI summary figure.

    volume: HU array (z, y, x); mask: boolean 3D ROI (z, y, x);
    data: 1-D array of ROI voxel HU values.
    """
    import matplotlib as mpl
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    z_sp, y_sp, x_sp = voxel_spacing
    n = len(data)
    mu, sd, med = float(np.mean(data)), float(np.std(data)), float(np.median(data))
    # Gaussian-based half-maximum interval: mean ± sqrt(2 ln 2) × SD
    lo_f, hi_f = mu - 1.177 * sd, mu + 1.177 * sd
    sel = data[(data >= lo_f) & (data <= hi_f)]
    fwhm_mean = float(np.mean(sel)) if len(sel) else 0.0

    # Display the sagittal and coronal planes that contain the most ROI voxels.
    x0 = int(np.argmax(mask.sum(axis=(0, 1))))
    y0 = int(np.argmax(mask.sum(axis=(0, 2))))
    views = [  # (title, image, ROI mask slice, horizontal mm/px, vertical mm/px)
        ("Sagittal", volume[::-1, :, x0], mask[::-1, :, x0], y_sp, z_sp),
        ("Coronal",  volume[::-1, y0, :], mask[::-1, y0, :], x_sp, z_sp)]

    rc = {'font.family': 'sans-serif',
          'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
          'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
          'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
          'axes.linewidth': 0.8, 'pdf.fonttype': 42, 'svg.fonttype': 'none'}
    with mpl.rc_context(rc):
        fig = Figure(figsize=(7.2, 3.9), dpi=dpi, facecolor='white')
        FigureCanvasAgg(fig)
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.8], left=0.02, right=0.985,
                              top=0.86, bottom=0.27, wspace=0.30)

        for k, (title, img, msk, h_sp, v_sp) in enumerate(views):
            ax = fig.add_subplot(gs[0, k])
            ax.imshow(img, cmap='gray', vmin=PUB_WINDOW[0], vmax=PUB_WINDOW[1],
                      aspect=v_sp / h_sp, interpolation='nearest')
            ax.contour(msk.astype(float), levels=[0.5], colors=[PUB_OUTLINE], linewidths=1.2)
            rows, cols = np.where(msk)
            pad_c, pad_r = 22 / h_sp, 22 / v_sp                 # 22 mm margin around the ROI
            ax.set_xlim(cols.min() - pad_c, cols.max() + pad_c)
            ax.set_ylim(rows.max() + pad_r, rows.min() - pad_r)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_anchor('N')          # align image panels at the top
            ax.set_title(title, pad=3)
            # 10 mm scale bar
            x_a, x_b = ax.get_xlim(); y_a, y_b = ax.get_ylim()
            bx = x_a + 0.06 * (x_b - x_a); by = y_a - 0.07 * (y_a - y_b)
            ax.plot([bx, bx + 10 / h_sp], [by, by], color='white', lw=2, solid_capstyle='butt')
            ax.text(bx + 5 / h_sp, by - 0.03 * (y_a - y_b), "10 mm", color='white',
                    ha='center', va='bottom', fontsize=8)
            if k == 0:
                ax.text(-0.02, 1.12, "A", transform=ax.transAxes, fontsize=12,
                        fontweight='bold', va='bottom', ha='left')

        axh = fig.add_subplot(gs[0, 2])
        w = np.full(n, 100.0 / n)
        axh.hist(data, bins=50, weights=w, color=PUB_BARS, alpha=0.75, edgecolor='white', linewidth=0.3)
        axh.set_ylim(0, axh.get_ylim()[1] * 1.55)     # head-room for the legend
        top = 0.66                                    # lines and shading stop below the legend
        axh.axvspan(lo_f, hi_f, ymax=top, color=PUB_SHADE, alpha=0.30, lw=0)
        axh.axvline(mu, ymax=top, color=PUB_MEAN, lw=1.6, ls='-')
        axh.axvline(med, ymax=top, color=PUB_MEDIAN, lw=1.6, ls='--')
        axh.axvline(fwhm_mean, ymax=top, color=PUB_FWHM_MEAN, lw=1.6, ls=':')
        axh.set_xlabel("HU"); axh.set_ylabel("Percentage of ROI voxels")
        axh.spines['top'].set_visible(False); axh.spines['right'].set_visible(False)
        axh.legend(handles=[
            Line2D([0], [0], color=PUB_MEAN, lw=1.6, ls='-', label=f"Mean HU ({mu:.1f})"),
            Line2D([0], [0], color=PUB_MEDIAN, lw=1.6, ls='--', label=f"Median HU ({med:.1f})"),
            Line2D([0], [0], color=PUB_FWHM_MEAN, lw=1.6, ls=':', label=f"FWHM mean HU ({fwhm_mean:.1f})"),
            Patch(facecolor=PUB_SHADE, alpha=0.30, label=f"Gaussian-based FWHM interval ({lo_f:.1f} to {hi_f:.1f})")],
            loc='upper right', frameon=False, handlelength=1.6)
        axh.text(-0.20, 1.12, "B", transform=axh.transAxes, fontsize=12,
                 fontweight='bold', va='bottom', ha='left')

        fig.text(0.5, 0.075,
                 f"3D ROI voxel count: {n:,}.  HU range: {float(np.min(data)):.0f} to "
                 f"{float(np.max(data)):.0f} (Gaussian-based FWHM interval: {lo_f:.1f} to {hi_f:.1f}).\n"
                 f"CT display window in panel A: {PUB_WINDOW[0]} to {PUB_WINDOW[1]} HU; "
                 f"orange line: 3D ROI outline.",
                 ha='center', va='center', fontsize=8, linespacing=1.5)
        fig.savefig(out_path, dpi=dpi, facecolor='white')
    return out_path


class WelcomeWindow:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Akdeniz Tool for HU Calculation")
        self.root.geometry("550x420")
        self.root.configure(bg="#002366")  # Akdeniz blue

        self.bg_canvas = tk.Canvas(self.root, width=550, height=420, 
                                   bg="#002366", highlightthickness=0)
        self.bg_canvas.pack(fill="both", expand=True)

        try:
            # Optional background logo, drawn at 30% opacity
            bg_img = Image.open(os.path.join(BASE_DIR, "background.png")).convert("RGBA")
            
            alpha = bg_img.split()[3]
            alpha = alpha.point(lambda p: p * 0.3) 
            bg_img.putalpha(alpha)
            
            bg_img = bg_img.resize((300, 230)) 
            
            self.bg_photo = ImageTk.PhotoImage(bg_img)
            self.bg_canvas.create_image(275, 210, image=self.bg_photo, anchor="center")
        except Exception:
            pass 

        title_font = ("Helvetica", 24, "bold")
        desc_font = ("Helvetica", 12)
        
        self.bg_canvas.create_text(275, 80, text="Akdeniz Tool for HU Calculation", 
                                   font=title_font, fill="#f37021", width=500, justify="center")

        short_desc = ("Hounsfield unit calculator using a 3D ROI\n\n"
                      "Calculates mean HU, median HU,\n"
                      "the Gaussian-based FWHM interval,\n"
                      "and FWHM mean HU.")
        
        self.bg_canvas.create_text(275, 180, text=short_desc, 
                                   font=desc_font, fill="#ecf0f1", justify="center")

        style = ttk.Style()
        style.theme_use('clam')
        style.configure("Akdeniz.TButton", font=("Helvetica", 12, "bold"), 
                        background="#f37021", foreground="white", bordercolor="#f37021")
        style.map("Akdeniz.TButton", background=[('active', '#d35400')])

        self.open_btn = ttk.Button(self.root, text="Open Spinal DICOM Folder", 
                                   style="Akdeniz.TButton", command=self.start_viewer)
        self.bg_canvas.create_window(275, 280, window=self.open_btn)

        footer_text = "Developed by Muhammad Ateya, Department of Neurosurgery, Akdeniz University"
        self.bg_canvas.create_text(275, 380, text=footer_text, 
                                   font=("Helvetica", 9, "italic"), fill="#bdc3c7")

        self.root.mainloop()

    def start_viewer(self):
        folder_path = filedialog.askdirectory(title="Select DICOM Folder")
        if folder_path:
            self.root.withdraw()
            VertebraHUViewer(folder_path, self.root)

class VertebraHUViewer:
    LABELS = [
        "LEFT: Posterior Superior", "LEFT: Anterior Superior",
        "LEFT: Anterior Inferior", "LEFT: Posterior Inferior",
        "RIGHT: Posterior Superior", "RIGHT: Anterior Superior",
        "RIGHT: Anterior Inferior", "RIGHT: Posterior Inferior"
    ]

    # 3D density-cloud display parameters (visualization only; they do not
    # affect ROI statistics). Larger TARGET_VOXELS / MAX_PTS_PER_GRADE give a
    # denser cloud but slower redraws.
    MIN_HU_DISPLAY = 100        # lower HU bound of the cloud (includes trabecular bone)
    HU_GRADES = 12              # number of HU bins between MIN_HU_DISPLAY and the local maximum
    TARGET_VOXELS = 220000      # approximate voxel count after regular subsampling
    MAX_PTS_PER_GRADE = 4500    # maximum plotted points per HU bin
    PREVIEW_PAD = 20            # voxels of padding around the ROI in the 3D views

    def __init__(self, folder_path, root_context):
        self.root = root_context
        self.labels = self.LABELS
        self.fig = None
        self._switching_subject = False
        if not self.load_subject(folder_path):
            self.root.deiconify()

    def load_subject(self, folder_path):
        """Load a DICOM folder and (re)build the viewer. Returns False on failure."""
        volume = self.load_dicom_series(folder_path)
        if volume is None:
            messagebox.showerror("Error", "No valid DICOM data found in that folder.")
            return False

        self.folder_path = folder_path
        self.volume = volume
        if not hasattr(self, 'voxel_spacing'):
            self.voxel_spacing = (1.0, 1.0, 1.0)
        self.points = []
        self.scatters, self.lines, self.patches = [], [], []
        self.stats_text = None
        self.zoom_mode = False
        self.show_3d_preview = False
        self.ax_3d_preview = None

        self.show_viewer()
        return True

    def load_dicom_series(self, folder_path):
        slices = []
        for filename in os.listdir(folder_path):
            filepath = os.path.join(folder_path, filename)
            try:
                ds = pydicom.dcmread(filepath, force=True)
                if not hasattr(ds, "PixelData"): continue
                img = ds.pixel_array.astype(np.float32)
                # Convert stored pixel values to HU: HU = RescaleSlope × pixel + RescaleIntercept
                img = img * float(getattr(ds, "RescaleSlope", 1)) + float(getattr(ds, "RescaleIntercept", 0))
                z_pos = float(ds.ImagePositionPatient[2]) if hasattr(ds, "ImagePositionPatient") else 0
                row_sp, col_sp = (float(v) for v in getattr(ds, "PixelSpacing", (1.0, 1.0)))
                slices.append((z_pos, img, row_sp, col_sp))
            except Exception:
                continue

        if not slices: return None
        slices.sort(key=lambda s: s[0])

        # Slice spacing (mm) is the median difference between consecutive
        # ImagePositionPatient z-positions; SliceThickness can differ from the
        # reconstruction increment.
        z_positions = [s[0] for s in slices]
        if len(z_positions) > 1:
            diffs = np.diff(z_positions)
            diffs = diffs[diffs > 0]
            z_spacing = float(np.median(diffs)) if len(diffs) else 1.0
        else:
            z_spacing = 1.0

        row_sp, col_sp = slices[0][2], slices[0][3]
        if any(abs(s[2] - row_sp) > 1e-3 or abs(s[3] - col_sp) > 1e-3 for s in slices):
            print("Warning: PixelSpacing is inconsistent across slices in this folder "
                  "(it may contain more than one series) - display proportions may be off.")

        # Voxel spacing (z, y, x) in mm, used for display proportions.
        self.voxel_spacing = (z_spacing, row_sp, col_sp)

        return np.stack([s[1] for s in slices])

    def new_subject(self):
        """Load a different DICOM folder into the running application."""
        folder_path = filedialog.askdirectory(title="Select DICOM Folder for New Subject")
        if not folder_path:
            self._focus_canvas()
            return

        old_fig = self.fig
        self._switching_subject = True
        loaded = self.load_subject(folder_path)
        if old_fig is not None:
            plt.close(old_fig)  # close while _switching_subject is True so the welcome window is not restored
        self._switching_subject = False
        if not loaded:
            # Loading failed: keep the current display.
            self._focus_canvas()

    def show_viewer(self):
        z_dim, y_dim, x_dim = self.volume.shape
        self.current_z, self.current_y, self.current_x = z_dim // 2, y_dim // 2, x_dim // 2

        self.fig, axes = plt.subplots(1, 3, figsize=(18, 8), facecolor='#002366')
        self.ax_ax, self.ax_cor, self.ax_sag = axes
        self.fig.canvas.manager.set_window_title('Akdeniz Tool for HU Calculation')

        for ax in axes:
            ax.set_facecolor('black')
            ax.tick_params(colors='white')

        # Display aspect ratios from voxel spacing: coronal and sagittal views
        # combine the slice axis (z) with an in-plane axis, so a unit aspect
        # would distort scans with non-cubic voxels.
        z_sp, y_sp, x_sp = self.voxel_spacing
        axial_aspect = y_sp / x_sp
        cor_aspect = z_sp / x_sp
        sag_aspect = z_sp / y_sp

        self.im_ax = self.ax_ax.imshow(self.volume[self.current_z], cmap="gray", aspect=axial_aspect)
        self.im_cor = self.ax_cor.imshow(self.volume[::-1, self.current_y, :], cmap="gray", aspect=cor_aspect)
        self.im_sag = self.ax_sag.imshow(self.volume[::-1, :, self.current_x], cmap="gray", aspect=sag_aspect)

        line_color = '#f37021'
        self.line_ax_h, = self.ax_ax.plot([0, x_dim], [self.current_y, self.current_y], color=line_color, lw=1.0, alpha=0.6, zorder=20)
        self.line_ax_v, = self.ax_ax.plot([self.current_x, self.current_x], [0, y_dim], color=line_color, lw=1.0, alpha=0.6, zorder=20)
        z_scr = z_dim - 1 - self.current_z
        self.line_cor_h, = self.ax_cor.plot([0, x_dim], [z_scr, z_scr], color=line_color, lw=1.0, alpha=0.6, zorder=20)
        self.line_cor_v, = self.ax_cor.plot([self.current_x, self.current_x], [0, z_dim], color=line_color, lw=1.0, alpha=0.6, zorder=20)
        self.line_sag_h, = self.ax_sag.plot([0, y_dim], [z_scr, z_scr], color=line_color, lw=1.0, alpha=0.6, zorder=20)
        self.line_sag_v, = self.ax_sag.plot([self.current_y, self.current_y], [0, z_dim], color=line_color, lw=1.0, alpha=0.6, zorder=20)

        self.stats_text = self.fig.text(0.02, 0.96,
                                        f"Step 1: {self.labels[0]}\n"
                                        "Left-Click: Move Crosshair | Right-Click: Place Point | "
                                        "Backspace: Undo | C: Clear | V: 3D Preview | Z: Zoom | N: New Subject",
                                        va='top', color='#f37021', fontsize=10, fontweight='bold',
                                        bbox=dict(facecolor='#002366', alpha=0.85, edgecolor='#f37021'))

        # Inset 3D ROI preview, hidden until toggled with V
        self.ax_3d_preview = self.fig.add_axes([0.40, 0.02, 0.20, 0.30], projection='3d')
        self.ax_3d_preview.set_facecolor('black')
        self.ax_3d_preview.set_visible(False)

        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.fig.canvas.mpl_connect('close_event', self.on_window_close)
        self.update_all()
        self._focus_canvas()
        plt.show()

    def _focus_canvas(self):
        """Give keyboard focus to the plot canvas."""
        try:
            self.fig.canvas.get_tk_widget().focus_force()
        except Exception:
            pass

    def on_window_close(self, event):
        # Do not restore the welcome window while switching subjects.
        if self._switching_subject:
            return
        self.root.deiconify()

    def update_all(self):
        z_dim, y_dim, x_dim = self.volume.shape
        self.im_ax.set_data(self.volume[self.current_z])
        self.im_cor.set_data(self.volume[::-1, self.current_y, :])
        self.im_sag.set_data(self.volume[::-1, :, self.current_x])

        z_scr = z_dim - 1 - self.current_z
        self.line_ax_h.set_ydata([self.current_y, self.current_y])
        self.line_ax_v.set_xdata([self.current_x, self.current_x])
        self.line_cor_h.set_ydata([z_scr, z_scr])
        self.line_cor_v.set_xdata([self.current_x, self.current_x])
        self.line_sag_h.set_ydata([z_scr, z_scr])
        self.line_sag_v.set_xdata([self.current_y, self.current_y])

        for obj in self.scatters + self.lines + self.patches:
            obj.remove()
        self.scatters, self.lines, self.patches = [], [], []

        for i, (py, pz_s, px) in enumerate(self.points):
            color = '#f37021' if i < 4 else '#00FF00'
            if px == self.current_x:
                s = self.ax_sag.scatter(py, pz_s, c=color, s=30, edgecolors='white', zorder=10)
                self.scatters.append(s)

        if len(self.points) == 8:
            pts = np.array(self.points)
            for start, color in [(0, '#f37021'), (4, '#00FF00')]:
                quad = Polygon(pts[start:start+4, 0:2], closed=True, color=color, alpha=0.25)
                self.ax_sag.add_patch(quad); self.patches.append(quad)
        
        controls_hint = ("Left-Click: Move Crosshair | Right-Click: Place Point | "
                          "Backspace: Undo | C: Clear | V: 3D Preview | Z: Zoom | N: New Subject")

        if self.zoom_mode:
            current_label = "Zoom Mode Active (Press Z to Resume Placing Points)"
            self.stats_text.set_text(current_label)
            self.stats_text.set_color('yellow')
        else:
            step_idx = len(self.points)
            if step_idx < 8:
                self.stats_text.set_text(f"Step {step_idx+1}: {self.labels[step_idx]}\n{controls_hint}")
            else:
                self.stats_text.set_text(f"Analysis Ready. Press Enter to Analyze or P for the Publication Figure\n{controls_hint}")
            self.stats_text.set_color('#f37021')

        if self.show_3d_preview:
            self.update_3d_preview()

        self.fig.canvas.draw_idle()

    def on_scroll(self, event):
        delta = 1 if event.button == 'up' else -1
        if event.inaxes == self.ax_ax: self.current_z = np.clip(self.current_z + delta, 0, self.volume.shape[0]-1)
        elif event.inaxes == self.ax_cor: self.current_y = np.clip(self.current_y + delta, 0, self.volume.shape[1]-1)
        elif event.inaxes == self.ax_sag: self.current_x = np.clip(self.current_x + delta, 0, self.volume.shape[2]-1)
        self.update_all()

    def on_click(self, event):
        self._focus_canvas()
        if self.zoom_mode or event.xdata is None or event.ydata is None:
            return

        z_dim, y_dim, x_dim = self.volume.shape

        if event.button == 3:  # Right-click: place an ROI point (sagittal view only)
            if event.inaxes != self.ax_sag:
                return
            if len(self.points) < 8:
                self.points.append((event.xdata, event.ydata, self.current_x))
                self.update_all()

        elif event.button == 1:  # Left-click: move crosshair
            if event.inaxes == self.ax_ax:
                self.current_x = int(np.clip(round(event.xdata), 0, x_dim - 1))
                self.current_y = int(np.clip(round(event.ydata), 0, y_dim - 1))
            elif event.inaxes == self.ax_cor:
                self.current_x = int(np.clip(round(event.xdata), 0, x_dim - 1))
                self.current_z = int(np.clip(z_dim - 1 - round(event.ydata), 0, z_dim - 1))
            elif event.inaxes == self.ax_sag:
                self.current_y = int(np.clip(round(event.xdata), 0, y_dim - 1))
                self.current_z = int(np.clip(z_dim - 1 - round(event.ydata), 0, z_dim - 1))
            else:
                return
            self.update_all()

    def on_key(self, event):
        if event.key == 'c':
            self.points = []
            self.update_all()
        elif event.key in ('backspace', 'u'):
            if self.points:
                self.points.pop()
                self.update_all()
        elif event.key == 'v':
            self.show_3d_preview = not self.show_3d_preview
            self.ax_3d_preview.set_visible(self.show_3d_preview)
            self.update_3d_preview()
            self.fig.canvas.draw_idle()
        elif event.key == 'z':
            self.zoom_mode = not self.zoom_mode
            if self.zoom_mode:
                self.fig.canvas.manager.toolbar.zoom()
            else:
                if self.fig.canvas.manager.toolbar.mode == 'zoom rect':
                    self.fig.canvas.manager.toolbar.zoom()
                elif self.fig.canvas.manager.toolbar.mode == 'pan/zoom':
                    self.fig.canvas.manager.toolbar.pan()
            self.update_all()
            self._focus_canvas()  # toolbar buttons can take keyboard focus
        elif event.key == 'enter' and len(self.points) == 8:
            self.analyze_plateau()
        elif event.key == 'p' and len(self.points) == 8:
            self.export_publication_figure()
        elif event.key == 'n':
            self.new_subject()

    def draw_density_cloud(self, ax, x_lo, x_hi, y_lo, y_hi, z_lo, z_hi):
        """Draw a graded HU point cloud for a voxel bounding box (visualization only).

        The crop is regularly subsampled and each HU bin is randomly thinned for
        display speed. This does not affect the ROI statistics, which use all ROI voxels.
        """
        sub = self.volume[z_lo:z_hi, y_lo:y_hi, x_lo:x_hi]

        # Subsampling step scales with crop size (display only).
        step = max(1, int(round((sub.size / self.TARGET_VOXELS) ** (1 / 3))))
        sub_s = sub[::step, ::step, ::step]

        local_max = float(sub_s.max()) if sub_s.size else self.MIN_HU_DISPLAY + 1
        if local_max <= self.MIN_HU_DISPLAY:
            local_max = self.MIN_HU_DISPLAY + 1
        edges = np.linspace(self.MIN_HU_DISPLAY, local_max, self.HU_GRADES + 1)
        cmap = plt.get_cmap('inferno')

        for i in range(self.HU_GRADES):
            lo, hi = edges[i], edges[i + 1]
            mask = (sub_s >= lo) & (sub_s <= hi) if i == self.HU_GRADES - 1 else (sub_s >= lo) & (sub_s < hi)
            zz, yy, xx = np.where(mask)
            if len(zz) == 0:
                continue
            if len(zz) > self.MAX_PTS_PER_GRADE:
                # Random thinning of this HU bin (display only)
                idx = np.random.choice(len(zz), self.MAX_PTS_PER_GRADE, replace=False)
                zz, yy, xx = zz[idx], yy[idx], xx[idx]
            frac = i / (self.HU_GRADES - 1)
            ax.scatter(xx * step + x_lo, yy * step + y_lo, zz * step + z_lo,
                       color=cmap(frac), alpha=0.05 + 0.6 * frac, s=3 + 8 * frac, linewidths=0)

    def update_3d_preview(self):
        if not self.show_3d_preview:
            return

        ax = self.ax_3d_preview
        ax.clear()
        ax.set_facecolor('black')
        ax.set_title("3D ROI Preview", color='white', fontsize=8)
        ax.tick_params(colors='white', labelsize=6)
        ax.grid(False)

        z_dim, y_dim, x_dim = self.volume.shape
        pad = self.PREVIEW_PAD

        if self.points:
            pts_arr = np.array(self.points, dtype=float)
            pts_v = pts_arr.copy()
            pts_v[:, 1] = z_dim - 1 - pts_v[:, 1]  # displayed (flipped) z -> volume z
            x_lo = int(max(0, pts_v[:, 2].min() - pad)); x_hi = int(min(x_dim, pts_v[:, 2].max() + pad))
            y_lo = int(max(0, pts_v[:, 0].min() - pad)); y_hi = int(min(y_dim, pts_v[:, 0].max() + pad))
            z_lo = int(max(0, pts_v[:, 1].min() - pad)); z_hi = int(min(z_dim, pts_v[:, 1].max() + pad))
        else:
            half = 25
            x_lo, x_hi = max(0, self.current_x - half), min(x_dim, self.current_x + half)
            y_lo, y_hi = max(0, self.current_y - half), min(y_dim, self.current_y + half)
            z_lo, z_hi = max(0, self.current_z - half), min(z_dim, self.current_z + half)

        # Ensure non-empty crop bounds
        if x_hi <= x_lo: x_hi = x_lo + 1
        if y_hi <= y_lo: y_hi = y_lo + 1
        if z_hi <= z_lo: z_hi = z_lo + 1

        self.draw_density_cloud(ax, x_lo, x_hi, y_lo, y_hi, z_lo, z_hi)

        if self.points:
            colors_pts = ['#f37021' if i < 4 else '#00FF00' for i in range(len(pts_v))]
            ax.scatter(pts_v[:, 2], pts_v[:, 0], pts_v[:, 1], c=colors_pts, s=25, edgecolors='white', zorder=10)

            if len(pts_v) == 8:
                for i in range(4):
                    j = (i + 1) % 4
                    ax.plot(pts_v[[i, j], 2], pts_v[[i, j], 0], pts_v[[i, j], 1], color='#f37021', lw=1.5)
                for i in range(4, 8):
                    j = 4 + (i + 1 - 4) % 4
                    ax.plot(pts_v[[i, j], 2], pts_v[[i, j], 0], pts_v[[i, j], 1], color='#00FF00', lw=1.5)
                for i in range(4):
                    ax.plot(pts_v[[i, i + 4], 2], pts_v[[i, i + 4], 0], pts_v[[i, i + 4], 1],
                           color='white', alpha=0.5, lw=1)

        # Box aspect from physical voxel spacing (mm), not voxel counts.
        z_sp, y_sp, x_sp = self.voxel_spacing
        ax.set_box_aspect(((x_hi - x_lo) * x_sp, (y_hi - y_lo) * y_sp, (z_hi - z_lo) * z_sp))
        ax.set_xlim(x_lo, x_hi); ax.set_ylim(y_lo, y_hi); ax.set_zlim(z_lo, z_hi)


    def roi_mask_and_data(self):
        """Build the 3D ROI mask and return (mask, ROI voxel HU values).

        Uses the same construction as analyze_plateau: for each sagittal slice x,
        the quadrilateral vertices are linearly interpolated between the first
        four points (first slice) and the last four points (second slice).
        """
        z_dim, y_dim, x_dim = self.volume.shape
        pts = np.array(self.points)
        pts_v = pts.copy()
        pts_v[:, 1] = z_dim - 1 - pts[:, 1]
        x_min, x_max = int(min(pts_v[:, 2])), int(max(pts_v[:, 2]))
        # (y, z) coordinate grid, identical for all x-slices
        yy, zz = np.meshgrid(np.arange(y_dim), np.arange(z_dim))
        grid_points = np.column_stack((yy.ravel(), zz.ravel()))
        mask3d = np.zeros(self.volume.shape, dtype=bool)
        for x in range(x_min, x_max + 1):
            f = (x - x_min) / (x_max - x_min) if x_max != x_min else 0
            quad = [(pts_v[i,0] + f*(pts_v[i+4,0]-pts_v[i,0]), pts_v[i,1] + f*(pts_v[i+4,1]-pts_v[i,1])) for i in range(4)]
            mask3d[:, :, x] = Path(quad).contains_points(grid_points).reshape(z_dim, y_dim)
        return mask3d, self.volume[mask3d]

    def export_publication_figure(self):
        mask3d, data = self.roi_mask_and_data()
        if data.size == 0:
            messagebox.showerror("Error", "The 3D ROI contains no voxels.")
            return
        path = filedialog.asksaveasfilename(title="Save publication figure", defaultextension=".png",
                                            initialfile="ROI_summary_figure.png",
                                            filetypes=[("PNG Image (300 dpi)", "*.png")])
        if path:
            make_publication_figure(self.volume, mask3d, self.voxel_spacing, data, path)
            pdf = os.path.splitext(path)[0] + ".pdf"
            make_publication_figure(self.volume, mask3d, self.voxel_spacing, data, pdf)
            messagebox.showinfo("Saved", f"Saved:\n{path}\n{pdf}")
        self._focus_canvas()

    def analyze_plateau(self):
        z_dim, y_dim, x_dim = self.volume.shape
        pts = np.array(self.points)
        pts_v = pts.copy()
        pts_v[:, 1] = z_dim - 1 - pts[:, 1]
        x_min, x_max = int(min(pts_v[:, 2])), int(max(pts_v[:, 2]))

        # 3D ROI: for each sagittal slice x, interpolate the quadrilateral vertices
        # between the first and last four points. The (y, z) grid is shared by all slices.
        yy, zz = np.meshgrid(np.arange(y_dim), np.arange(z_dim))
        grid_points = np.column_stack((yy.ravel(), zz.ravel()))

        all_voxels = []
        for x in range(x_min, x_max + 1):
            f = (x - x_min) / (x_max - x_min) if x_max != x_min else 0
            quad = [(pts_v[i,0] + f*(pts_v[i+4,0]-pts_v[i,0]), pts_v[i,1] + f*(pts_v[i+4,1]-pts_v[i,1])) for i in range(4)]
            mask = Path(quad).contains_points(grid_points).reshape(z_dim, y_dim)
            all_voxels.extend(self.volume[mask, x])

        if not all_voxels: return
        data = np.array(all_voxels)
        voxel_count = len(data)
        
        mu, std = np.mean(data), np.std(data)
        med = np.median(data)
        vmin, vmax = np.min(data), np.max(data)
        
        # Gaussian-based half-maximum interval: mean ± sqrt(2 ln 2) × SD
        plateau_min = mu - 1.177 * std
        plateau_max = mu + 1.177 * std
        
        fwhm_mask = (data >= plateau_min) & (data <= plateau_max)
        fwhm_data = data[fwhm_mask]
        fwhm_mean = np.mean(fwhm_data) if len(fwhm_data) > 0 else 0

        # Clipboard output: tab-separated for pasting into a spreadsheet
        cb_header = "Voxel Count\tOverall Mean\tOverall Median\tFWHM Mean\tFWHM Min\tFWHM Max\tMin\tMax"
        cb_stats = f"{voxel_count}\t{mu:.2f}\t{med:.2f}\t{fwhm_mean:.2f}\t{plateau_min:.2f}\t{plateau_max:.2f}\t{vmin:.2f}\t{vmax:.2f}"
        results_clipboard = f"{cb_header}\n{cb_stats}"
        
        # On-screen output: space-aligned
        disp_header = f"{'Voxel Count':<14}{'Overall Mean':<16}{'Overall Median':<16}{'FWHM Mean':<14}{'FWHM Min':<14}{'FWHM Max':<14}{'Min':<10}{'Max':<10}"
        disp_stats = f"{str(voxel_count):<14}{mu:<16.2f}{med:<16.2f}{fwhm_mean:<14.2f}{plateau_min:<14.2f}{plateau_max:<14.2f}{vmin:<10.1f}{vmax:<10.1f}"
        results_display = f"{disp_header}\n{disp_stats}\n(Copied to Clipboard)"

        try:
            self.root.clipboard_clear()
            self.root.clipboard_append(results_clipboard)
            self.root.update()
        except Exception:
            pass

        fig = plt.figure(figsize=(18, 7), facecolor='#002366')
        fig.canvas.manager.set_window_title("Akdeniz Tool for HU Calculation - 3D ROI and Statistics")
        fig.canvas.mpl_connect('close_event', lambda evt: self._focus_canvas())
        
        ax1 = fig.add_subplot(131)
        ax1.set_facecolor('black')
        # 10-HU bins; gold bins overlap the Gaussian-based FWHM interval
        clusters = (data // 10) * 10
        unique_cl, counts = np.unique(clusters, return_counts=True)
        percentages = (counts / len(data)) * 100
        colors = ['#FFD700' if (c >= plateau_min - 10 and c <= plateau_max) else '#4682B4' for c in unique_cl]
        ax1.bar(unique_cl, percentages, width=8, color=colors, alpha=0.9, edgecolor='white', lw=0.5)
        ax1.set_title("HU Distribution by 10-HU Bin\n(Gold = Gaussian-Based FWHM Interval)", color='white', fontsize=10)
        ax1.set_xlabel("HU (10-HU Bins)", color='white')
        ax1.set_ylabel("Percentage of ROI Voxels", color='white')
        ax1.tick_params(colors='white')

        ax2 = fig.add_subplot(132)
        ax2.set_facecolor('black')
        ax2.hist(data, bins=50, density=True, alpha=0.2, color='#bdc3c7')
        x_axis = np.linspace(vmin, vmax, 200)
        y_gauss = norm.pdf(x_axis, mu, std)
        ax2.plot(x_axis, y_gauss, color='white', lw=1, alpha=0.5, ls='--')
        
        # Shade the Gaussian-based FWHM interval
        ax2.fill_between(x_axis, 0, y_gauss, where=(x_axis >= plateau_min) & (x_axis <= plateau_max), color='#FFD700', alpha=0.4)
        
        # Summary statistics
        ax2.axvline(mu, color='#FF4500', lw=2, label=f'Mean HU: {mu:.1f}')
        ax2.axvline(med, color='#00FF00', lw=1.5, ls='-.', label=f'Median HU: {med:.1f}')
        ax2.axvline(fwhm_mean, color='#00FFFF', lw=2, label=f'FWHM Mean HU: {fwhm_mean:.1f}')
        
        # Gaussian-based FWHM interval limits
        ax2.axvline(plateau_min, color='#FFD700', lw=1, ls='--', label=f'FWHM Interval Min: {plateau_min:.1f}')
        ax2.axvline(plateau_max, color='#FFD700', lw=1, ls='--', label=f'FWHM Interval Max: {plateau_max:.1f}')
        
        ax2.axvline(vmin, color='white', lw=1, ls=':', label=f'ROI Min: {vmin:.0f}')
        ax2.axvline(vmax, color='white', lw=1, ls=':', label=f'ROI Max: {vmax:.0f}')
        
        ax2.set_title("HU Distribution of the 3D ROI", color='white')
        ax2.legend(fontsize='7', loc='upper right', facecolor='#2c3e50', labelcolor='white')
        ax2.tick_params(colors='white')

        ax3 = fig.add_subplot(133, projection='3d')
        ax3.set_facecolor('black')
        pts_3d = pts_v.copy()

        pad = self.PREVIEW_PAD
        bx_lo = int(max(0, pts_3d[:, 2].min() - pad)); bx_hi = int(min(x_dim, pts_3d[:, 2].max() + pad))
        by_lo = int(max(0, pts_3d[:, 0].min() - pad)); by_hi = int(min(y_dim, pts_3d[:, 0].max() + pad))
        bz_lo = int(max(0, pts_3d[:, 1].min() - pad)); bz_hi = int(min(z_dim, pts_3d[:, 1].max() + pad))
        if bx_hi <= bx_lo: bx_hi = bx_lo + 1
        if by_hi <= by_lo: by_hi = by_lo + 1
        if bz_hi <= bz_lo: bz_hi = bz_lo + 1
        self.draw_density_cloud(ax3, bx_lo, bx_hi, by_lo, by_hi, bz_lo, bz_hi)

        for i in range(4):
            j = (i + 1) % 4
            ax3.plot(pts_3d[[i, j], 2], pts_3d[[i, j], 0], pts_3d[[i, j], 1], color='#f37021', lw=2)
        for i in range(4, 8):
            j = 4 + (i + 1 - 4) % 4
            ax3.plot(pts_3d[[i, j], 2], pts_3d[[i, j], 0], pts_3d[[i, j], 1], color='#00FF00', lw=2)
        for i in range(4):
            ax3.plot(pts_3d[[i, i+4], 2], pts_3d[[i, i+4], 0], pts_3d[[i, i+4], 1], color='white', alpha=0.6)
        ax3.set_title(f"3D ROI Shape\n({voxel_count:,} Voxels)", color='white')
        ax3.tick_params(colors='white', labelsize=7)
        ax3.xaxis.label.set_color('white'); ax3.yaxis.label.set_color('white'); ax3.zaxis.label.set_color('white')
        z_sp, y_sp, x_sp = self.voxel_spacing
        ax3.set_box_aspect(((bx_hi - bx_lo) * x_sp, (by_hi - by_lo) * y_sp, (bz_hi - bz_lo) * z_sp))
        ax3.grid(False)

        self.stats_text.set_text(results_display)
        self.stats_text.set_fontsize(8)
        self.stats_text.set_fontfamily('monospace')
        
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    WelcomeWindow()