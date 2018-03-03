import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from seaborn import xkcd_rgb

class ProteinChordDiagram:
    """
    Makes and saves a Circos-like 'chord diagram' for a polymer, where the lines represent monomer-monomer contacts, and a 'ChiP-Seq' profile is also plotted as a part of the same figure.
    """
    def __init__(self, sites, polymer_length, radius = 1, max_radius=0.9, color_scheme='dark', size=600):
        self.color_scheme = color_scheme
        self.attrs = group.attrs
        self.L = polymer_length
        self.radius = radius
        self.max_radius = max_radius
        self.set_colors()
        self.global_rotation = np.pi/2
        self.dpi = 150
        self.default_size = size
        self.sites = sites

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        plt.clf()
        plt.close()

    def set_colors(self):
        """ Sets the color scheme, based on the private variable self.color_scheme """
        assert self.color_scheme in ('light', 'dark')
        if self.color_scheme == 'dark':
            self.bgcolor = (.1, .1, .1)
            self.color_text = (1,1,1)
            self.color_polymer = (.8, .8, .8)
            self.color_protein_bivalent = xkcd_rgb['cerulean']
            self.color_protein_monovalent = xkcd_rgb['turquoise']
            self.color_ori = xkcd_rgb['green']
            self.color_ter = xkcd_rgb['light red']
        else:
            self.bgcolor = (1,1,1)
            self.color_text = (0,0,0)
            self.color_polymer = (.5, .5, .5)
            self.color_protein_bivalent = xkcd_rgb['darkish red']
            self.color_protein_monovalent = xkcd_rgb['turquoise']
            self.color_ori = xkcd_rgb['green']
            self.color_ter = xkcd_rgb['darkish red']

    def initialize(self):
        """ Sets the plot style defaults and creates a new figure and axes.  """
        plt.figure(facecolor=self.bgcolor, edgecolor='none', figsize=(self.default_size/self.dpi,self.default_size/self.dpi), dpi=self.dpi)
        plt.subplot(111, polar=True)
        plt.axis('off')

    def get_protein_sites(self):
        """
        Returns the *only* the protein sites_all (i.e. not their coordinates). Return
        data does not include the conformation of the polymer (i.e. whether it's
        a monovalent, bivalent protein etc.).
        :returns: Nx2 array of integers that specify the sites_all of the various proteins
        on the polymer. The two columns specify the sites_all.
        """
        return self.sites  # only up to the second column because the third column is the 'conformation' of the protein, which is information we don't need for the visualization

    def compute_angles(self, _sites):
        """
        Returns the angles at which certain sites_all along the polymer lie.
        :param _sites: Numpy array of integers in [0, self.L)
        :return: floats of the angles of the same size as :param sites_all
        """
        return np.mod(2 * np.pi * _sites / self.L + self.global_rotation, 2*np.pi)

    def compute_central_radius_of_curvature(self, angle_1, angle_2):
        """
        To smoothly connect a curve between two coordinates that is also orthogonal to the circle on which those two
        coordinates lie we use an intermediate supporting point. This function computes the 'radius'. The more straight
        the connecting line runs, the higher this radius, and the more curved the line is, the smaller this radius. If
        the connecting line runs completely straight (over an angle of pi), the connecting radius is defined as self.radius.
        :param angle_1: Angle of the one coordinate
        :param angle_2: Angle of the other coordinate
        :return: Angle of the supporting coordinate
        """
        deltaAngle = abs(angle_2 - angle_1) % (2 * np.pi)
        return self.radius * deltaAngle / np.pi

    def cart2pol(self, point):
        """
        Converts Cartesian coordinates to polar coordinates
        :param point: Tuple of x,y coordinates
        :return: Tuple of (angle, radius)
        """
        x, y = point[0], point[1]
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return phi + self.global_rotation, rho

    @staticmethod
    def compute_smallest_loop(site_1, site_2, max_size):
        """
        Returns the smallest connecting distance between site_1, site_2, given that the maximal distance is max_size
        :param site_1:
        :param site_2:
        :param max_size: Maximum connecting distance
        :return:
        """
        result = abs(site_2 - site_1)
        if result >= max_size/2:
            result = max_size - result
        return result

    def compute_central_bezier_point(self, angle_1, angle_2):
        """
        To smoothly connect a curve between two coordinates that is also orthogonal to the circle on which those two
        coordinates lie we use an intermediate supporting point. This function returns a tuple of x,y coordinates of
        this supporting point. These x,y coordinates lie on the Bezier curve.
        :param angle_1: Angle of the one coordinate
        :param angle_2: Angle of the other coordinate
        :return: x,y coordinates of the supporting coordinate
        """
        loop_size = self.compute_smallest_loop(angle_1, angle_2, 2*np.pi)
        centralRadius = (self.radius - loop_size/np.pi)/1.15

        centralAngle = (angle_1 + angle_2) / 2  # angle at which the supporting point is
        if not np.isclose(min(angle_1, angle_2) + loop_size, max(angle_1, angle_2)): # if the two smallest distance between the two angles overlaps with the origin region....
            centralAngle -= np.pi # then rotate it by pi (e.g. for angles (0+0.1,2*pi-0.1) the central angle would be pi, whereas it should be identically 0)
        centralRadius = self.max_radius if centralRadius > self.max_radius else centralRadius
        return centralAngle, centralRadius

    def draw_vertex_from_polar(self, point_1, point_2, lw=1, alpha=1.):
        """
        Draws a connecting line between two sites_all on the canvas based on the polar points point_1=(phi,rho), point_2=(phi,rho)
        :param point_1: Tuple or list of phi,ro coordinates
        :param point_2: Tuple or list of phi,ro coordinates
        :param lw: Line width of the vertex drawn on the canvas
        :param alpha: Opacity of lines
        """
        headAngle, _ = point_1
        tailAngle, _ = point_2
        centralPoint = self.compute_central_bezier_point(headAngle, tailAngle)  # central point on the Bezier curve
        verts = [
            point_1,  # P0
            centralPoint,  # P1
            centralPoint,  # P2
            point_2,  # P3
        ]
        codes = [Path.MOVETO,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 ]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', lw=lw, ec=self.color_protein_bivalent, alpha=alpha)
        plt.gca().add_patch(patch)

    def draw_ori(self, radius = 1, interpolation_ratio=10, alpha=1, ori_size=1):
        """
        Draws a red dot at site=0 to represent the origin region.
        """
        thermal_variance = ori_size

        if thermal_variance < (self.L / 2) ** 2:  # if this condition does not hold, there's aspecific binding, and it would be senseless to draw an origin
            delta = int(np.sqrt(thermal_variance))
            begin, end = int(-delta), int(delta+1)
            begin, end = 0, 2
            sites = np.linspace(begin, end, num=interpolation_ratio*(end-begin))
            plt.plot(self.compute_angles(sites), len(sites)*[radius], c=self.color_ori, alpha=alpha, lw=3, zorder=100)

    def draw_ter(self, radius = 1, alpha=1, interpolation_ratio=10, ter_size=0):
        """
        Draws a red dot at site=0 to represent the origin region.
        :param radius:
        :param alpha:
        :return:
        """
        thermal_variance = ter_size
        if thermal_variance < (self.L / 2) ** 2:  # if this condition does not hold, there's aspecific binding, and it would be senseless to draw an origin
            delta = int(np.sqrt(thermal_variance))
            begin, end = int(self.L/2-delta), int(self.L/2+delta+1)
            angles = np.linspace(begin, end, num=interpolation_ratio*(end-begin))
            plt.plot(self.compute_angles(angles), len(angles)*[radius], c=self.color_ter, alpha=alpha, lw=3, zorder=100)

    def draw_bounding_circle(self, lw=1):
        """
        Draws a thick line to represent the 'genome'.
        :param lw:
        :return:
        """
        angles = np.arange(0, 2*np.pi, step=0.01)
        plt.plot(angles, self.radius * np.ones(len(angles)), c=self.color_polymer, lw=lw, zorder=-999)

    def draw_vertices(self, lw=2, alpha=1):
        """
        Successively calls the method drawVertex to draw multiple connecting Bezier curves between sites_all.
        :param lw: Line-width of the Bezier curve
        :param alpha: opacity of lines
        """
        data = self.get_protein_sites()
        for sites in data:
            angle_1, angle_2 = self.compute_angles(sites)
            if sites[0] != sites[1]:
                self.draw_vertex_from_polar((angle_1, self.radius), (angle_2, self.radius), lw=lw, alpha=alpha)
                for angle in (angle_1, angle_2):

                    plt.plot(angle, self.radius, 'o', c=self.color_protein_bivalent, markersize=8*np.sqrt(160/(self.L+80)), alpha=alpha)
            else:
                plt.plot(angle_1, self.radius, 'o', c=self.color_protein_monovalent, markersize=8*np.sqrt(160/(self.L+80)), alpha=alpha)

    def add_text(self):
        """
        Adds axis labels to the current canvas
        """
        plt.annotate(r'$i=0$', xy=(np.pi/2, self.radius+0.08), color=self.color_text, horizontalalignment='center')
        plt.annotate(r'$i=\frac{1}{2}N_m$', xy=(-np.pi/2, self.radius+0.09), color=self.color_text, horizontalalignment='center', verticalalignment='top')

    def draw_all(self, add_text=True, annotate='', fig_init_handle=None, alpha=1, lw=2):
        """ Draw everything on the canvas; protein positions, vertices, etc.
        :param add_text: Boolean that specifies whether or not to add a phi-axis label.
        :param annotate: Extra text to add to the top of the canvas
        """
        if fig_init_handle is None:
            self.initialize()
        else:
            pass
        self.draw_bounding_circle()
        self.draw_vertices(lw=lw, alpha=alpha)
        self.draw_ori()
        self.draw_ter()
        if add_text:
            self.add_text()
        print('annotate text is:', annotate)
        plt.gca().annotate(annotate, xy=(np.pi / 2, 1.20), horizontalalignment='center', fontsize=6, color = (1,1,1) if self.color_scheme=='dark' else (0,0,0))
        plt.gca().set_rmax(1.50)  # make the radius a definite value so that view of the axes limits
        plt.tight_layout()

    def savefig(self, output_file):
        if self.color_scheme == 'dark':
            plt.savefig(output_file, facecolor=self.bgcolor, edgecolor='none', dpi=self.dpi)
        else:
            plt.savefig(output_file, edgecolor='none', dpi=self.dpi, transparent=True)

if __name__ == '__main__':
    import h5py
    file_name = '/Users/C.Miermans/projects/chromosome/python_packages/chromosome_visualization/test_data/raw_data.h5'
    File = h5py.File(file_name, mode='r')
    group = File[list(File)[0]]
    sites = np.array(group[f'{list(group)[-1]}/proteins'])[:, :2]

    X = ProteinChordDiagram(sites, group.attrs['n_monomers'][0], color_scheme='light')
    # X.draw_vertices(list(group)[0])
    X.draw_all(add_text=True, alpha=0.5)

    X.savefig('example.png')
    # X.draw_and_save_sequence(times, output_dir='../test_data/')
    # X.plot('chord_diagram')