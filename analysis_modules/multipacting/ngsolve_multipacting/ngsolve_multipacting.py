import matplotlib.pyplot as plt
import copy
import matplotlib
import numpy as np

plt.rcParams["figure.figsize"] = (9, 4)

q0 = 1.60217663e-19
m0 = 9.1093837e-31
mu0 = 4*np.pi*1e-7
eps0 = 8.85418782e-12
c0 = 299792458

import matplotlib
from matplotlib import cm

% matplotlib
notebook
# calling it a second time may prevent some graphics errors
% matplotlib
notebook
plt.close()
plt.rcParams["figure.figsize"] = (9, 4)
import copy
import matplotlib
from matplotlib import cm


def get_cmap(n, name='jet'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return cm.get_cmap(name, n)


#     return matplotlib.colormaps[name]

def cross(a, b):
    c1 = np.array(a)[:, 1].real * np.array(b)[:, 0].real
    c2 = -np.array(a)[:, 0].real * np.array(b)[:, 0].real
    return np.array([c1, c2]).T


def segment_intersection(line1, line2):
    x1, y1 = line1[0]
    x2, y2 = line1[1]
    x3, y3 = line2[0]
    x4, y4 = line2[1]

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if denom == 0:
        return False, (0, 0)

    tt = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    uu = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / denom
    if (0 <= tt <= 1) and (0 <= uu <= 1):
        px, py = x1 + tt * (x2 - x1), y1 + tt * (y2 - y1)
        return True, (px, py)
    else:
        return False, (0, 0)


def get_neighbours(surf_pts, idx):
    surf_pts_neigs = np.array(surf_pts[idx])
    #     print("indx", surf_pts_in_box_and_neigs)
    return surf_pts_neigs[surf_pts_neigs[:, 0].argsort()]


class Particles:
    def __init__(self, xrange, init_v, bounds):
        self.bounds = np.array(bounds)
        self.x = self.bounds[(self.bounds[:, 0] > xrange[0]) & (self.bounds[:, 0] < xrange[1])]

        # get normal pointing inwards of emission points
        self.pt_normals = np.ones(self.x.shape)
        # get neighbouring points
        res, idxs = self.distance(1)
        for nn, idx in enumerate(idxs):
            n12 = self.get_point_normal(idx)
            self.pt_normals[nn] = n12 / np.linalg.norm(n12)
        # convert velocity from eV to m/s
        init_v = np.sqrt(2 * init_v * q0 / m0)
        self.u = self.pt_normals * init_v

        self.dist = np.zeros(self.x.shape)
        self.collision_state = np.zeros(self.x.shape)
        cmap = get_cmap(len(self.x))
        self.colors = np.array([cmap(i) for i in range(len(self.x))])
        self.len = len(self.x)

        self.x_old, self.u_old = np.zeros(self.x.shape), np.zeros(self.u.shape)  # <- hold solution from prev time step
        self.x_temp, self.u_temp = np.zeros(self.x.shape), np.zeros(
            self.u.shape)  # <- hold tentative change in position

        self.record = [self.x]
        self.lost_particles = []

    def update(self, e, b):
        self.x_old = copy.deepcopy(self.x)
        self.u_old = copy.deepcopy(self.u)

        self.x_temp = copy.deepcopy(self.x)
        self.u_temp = copy.deepcopy(self.u)

        self.u = self.u + q0 / m0 * (e + cross(self.u, b)) * dt
        self.x = self.x + self.u * dt

    def distance(self, n):
        closest_dist_to_surface = []
        indx = []
        for xx in self.x:
            res = [np.linalg.norm(xx - xsrf) for xsrf in self.bounds]
            closest_dist_to_surface.append(np.min(res))
            indx.append([list(np.argpartition(res, n)[0:n])][0])  # <- index of two closest surface points

        return np.atleast_2d(closest_dist_to_surface).T, indx

    def remove(self, ind):
        self.x = np.delete(self.x, ind, axis=0)
        self.u = np.delete(self.u, ind, axis=0)

        self.x_old = np.delete(self.x_old, ind, axis=0)
        self.u_old = np.delete(self.u_old, ind, axis=0)

        self.x_temp = np.delete(self.x_temp, ind, axis=0)
        self.u_temp = np.delete(self.u_temp, ind, axis=0)

        self.colors = np.delete(self.colors, ind, axis=0)

        self.len = len(self.x)

    def colors(self):
        return self.colors

    def get_point_normal(self, idx):
        # calculate normal as average of connecting edget normals
        # assumpution is that the surface points are ordered in increasing x
        x0, y0 = self.bounds[idx[0] - 1]
        x1, y1 = self.bounds[idx[0]]
        x2, y2 = self.bounds[idx[0] + 1]

        dx1, dy1 = x1 - x0, y1 - y0
        dx2, dy2 = x2 - x1, y2 - y1

        n1 = -np.array([-dy1, dx1])
        n2 = -np.array([-dy2, dx2])
        n12 = n1 + n2

        return n12

    def update_record(self):
        self.record.append(self.x)


# create plot
fig, ax = plt.subplots()

n = 3
q0 = 1.60217663e-19
m0 = 9.1093837e-31

# get surface points
pec_boundary = mesh.Boundaries("default")
bel = [xx.vertices for xx in pec_boundary.Elements()]
bel_unique = list(set(itertools.chain(*bel)))
xpnts_surf = sorted([mesh.vertices[xy.nr].point for xy in bel_unique])
xsurf = np.array(xpnts_surf)
# print(np.array(xpnts_surf))

# # check all surface points
# checkk = np.array([np.array([gfu_E[n](mesh(*p))]) for p in xsurf])
# print('check', checkk)

# get xmin, xmax, ymin

xmin = face.vertices.Min(X).p[0]
xmax = face.vertices.Max(X).p[0]
ymin = 0

# define one internal point
xin = [0, 0]

# xpnts_surf = cav_geom[(cav_geom[0] > 0) & (cav_geom[1] > min(cav_geom[1])) & (cav_geom[1] < max(cav_geom[1]))]
# xsurf = xpnts_surf.to_numpy()

ax.plot(particles.bounds[:, 0], particles.bounds[:, 1], marker='o', ms=3, mec='k', mfc='none')

no_of_remaining_particles = []
Epks_array = np.linspace(1, 91, 31)
# Epks_array = [52]
for epk in Epks_array:
    scl = 1 / Epk * 1e6 * epk  # <- scale Epk to 1 MV/m and multiply by sweep value

    v_init = 200  # eV
    # v_init = [(), ()]  # <- initial velocities in eV
    particles = Particles([-0.05, 0.05], v_init, xsurf)
    print("No of initial particles: ", particles.len)
    print("= " * 50)

    # gfu_E = sclgfu_E[n])
    # gfu_H = scl*gfu_H[n]

    t = 0
    dt = 1e-11
    w = 2 * np.pi * freq_fes[n] * 1e6
    PHI = [0]  # <- initial phase

    # move particles with initial velocity. ensure all initial positions after first move lie inside the bounds
    particles.x = particles.x + particles.u * dt
    e_check_list = []
    time_list = []
    error = False
    counter = 0

    # calculate time for 10 cycles, 20 alternations
    T = 1 / (freq_fes[n] * 1e6) * 10
    #     print("Time for 10T", T)

    PLOT = False
    # sweep Epk

    for phi in PHI:
        while t < T and particles.len != 0:
            #         if counter%10 == 0:
            #             PLOT = True
            #         else:
            #             PLOT = False

            # get e and b field from eigenmode analysis at particle current position
            e = []
            for pind, p in enumerate(particles.x):
                try:
                    e.append(gfu_E[n](mesh(*p)))
                except Exception as exc:
                    error = True
                    print(pind, particles.x[pind])
                    #                 print(pind+1, particles.x[pind+1])
                    #                 print(pind+2, particles.x[pind+2])
                    #                 print(pind+3, particles.x[pind+3])
                    #                 print(pind+4, particles.x[pind+4])
                    print(exc)
                    break
            if error:
                break

            e = scl * np.array(e) * exp(1j * (w * t + phi))
            #         e = scl*np.array([gfu_E[n](mesh(*p)) for p in particles.x])*exp(1j*(w*t + phi))
            b = mu0 * scl * np.array([gfu_H[n](mesh(*p)) for p in particles.x]) * exp(1j * (w * t + phi))

            #         # check field at single point to see when it changes to zero
            #         e_check = scl*np.array(gfu_E[n](mesh(*(0, 0.1))))*exp(1j*(w*t + phi))
            # #         print(np.linalg.norm(e_check.real))
            #         e_check_list.append(e_check.real)
            #         time_list.append(t)

            particles.update(e.real, b.real)
            ############################################################################
            for col_ind, ptts in enumerate(particles.x):
                #             plt.scatter(ptts[0], ptts[1], color=particles.colors[col_ind], s=5)
                # calculate velocity magnitude
                vmag = np.linalg.norm(particles.u[col_ind])
                #             print(vmag)
                if PLOT:
                    ax.scatter(ptts[0], ptts[1], color=particles.colors[col_ind], s=5)
            #############################################################################

            # check if particle close to boundary
            res, indx = particles.distance(100)
            #         print(res, len(res))
            lost_particles_indx = []
            for ind, (r, idx) in enumerate(zip(res, indx)):
                if r < 5e-2:  # point at boundary, calculate new field value
                    #                 print("res les", len(res), particles.len, ind)
                    # check if point is inside or outside of region
                    # get surface points neighbours
                    surf_pts_neigs = get_neighbours(xsurf, idx)
                    #                 print('aftersort', surf_pts_neigs)

                    # check for intersection
                    #                 dp = False  # <- flag to check if a point has been deleted
                    for pt_neigs_ind in range(1, len(surf_pts_neigs)):
                        line1 = (particles.x[ind], xin)
                        line2 = (surf_pts_neigs[pt_neigs_ind - 1], surf_pts_neigs[pt_neigs_ind])

                        bool_intc, (x_intc, y_intc) = segment_intersection(line1, line2)
                        #                     print('\t ', pt_neigs_ind, bool_intc)

                        if bool_intc:

                            # plot line to
                            if PLOT:
                                ax.plot([particles.x[ind][0], xin[0]], [particles.x[ind][1], xin[1]], c='cyan')
                                ax.plot([surf_pts_neigs[pt_neigs_ind - 1][0], surf_pts_neigs[pt_neigs_ind][0]],
                                        [surf_pts_neigs[pt_neigs_ind - 1][1], surf_pts_neigs[pt_neigs_ind][1]],
                                        color='r')
                                #                             print([surf_pts_neigs[pt_neigs_ind-1][0], surf_pts_neigs[pt_neigs_ind][0]],
                                #                                      [surf_pts_neigs[pt_neigs_ind-1][1], surf_pts_neigs[pt_neigs_ind][1]])
                                #                             plot interseciton point
                                #                             print("plotting intersection", x_intc, y_intc)

                                ax.scatter(x_intc, y_intc, color='k', s=100, marker='P')

                            ####################################################
                            # get intersection with old point. loop through points again.
                            # the surface edge a line between an outside point and the origin intersects
                            # might be different from that with which the line between old and new point intersects

                            line11 = (particles.x[ind], particles.x_old[ind])
                            for pt_neigs_ind2 in range(1, len(surf_pts_neigs)):
                                line22 = (surf_pts_neigs[pt_neigs_ind2 - 1], surf_pts_neigs[pt_neigs_ind2])
                                bool_intc_p, (x_intc_p, y_intc_p) = segment_intersection(line11, line22)
                                if bool_intc_p:
                                    break

                            #                         print('bool_intc_p', bool_intc_p, line11, line22)
                            if bool_intc_p:
                                if PLOT:
                                    ax.plot([particles.x[ind][0], particles.x_old[ind][0]],
                                            [particles.x[ind][1], particles.x_old[ind][1]], c='b')

                                    ax.plot([surf_pts_neigs[pt_neigs_ind2 - 1][0], surf_pts_neigs[pt_neigs_ind2][0]],
                                            [surf_pts_neigs[pt_neigs_ind2 - 1][1], surf_pts_neigs[pt_neigs_ind2][1]],
                                            color='g')

                                    ax.scatter(x_intc_p, y_intc_p, color='y', s=100, marker='P')
                            ####################################################

                            ###################################################################
                            # calculate time fraction required to hit wall and calculate the field values and
                            # velocity at this time fraction
                            # calculate parameter tfrac
                            dt_frac = np.linalg.norm([x_intc_p, y_intc_p] - particles.x_old[ind]) / np.linalg.norm(
                                particles.x[ind] - particles.x_old[ind])
                            t_frac = t - dt * (1 - dt_frac)
                            #                         print(dt_frac, t_frac, t)

                            #  calculate field values at this time which is a (fraction of dt) + t
                            e = scl * np.array([gfu_E[n](mesh(*(x_intc_p, y_intc_p)))]) * exp(1j * (w * t_frac + phi))
                            b = mu0 * scl * np.array([gfu_H[n](mesh(*(x_intc_p, y_intc_p)))]) * exp(
                                1j * (w * t_frac + phi))

                            #                         particles.u_temp[ind] = particles.u[ind] + q0/m0*(e.real + cross([particles.u[ind]], b.real))*dt*dt_frac
                            # check if the e-field surface normal is close to zero indicating a possible change in field
                            particles.x_temp[ind] = particles.x_old[ind] + particles.u[ind] * dt * dt_frac
                            #                         print(particles.x_temp[ind], particles.x_old[ind], particles.x[ind])

                            # x_intc, y_intc is found on line22. get line normal.
                            # sort line to make sure x is increasing
                            line22 = np.array(line22)
                            line22 = line22[line22[:, 0].argsort()]
                            line22_normal = -np.array([-(line22[1][1] - line22[0][1]), line22[1][0] - line22[0][0]])
                            #                         line22_normal = line22_normal/np.linalg.norm(line22_normal)
                            #                         if PLOT:
                            #                             # plot surface norm
                            #                             ax.quiver(*[x_intc_p, y_intc_p], *line22_normal, scale=15)
                            #                             # plot efield
                            #                             ax.quiver(x_intc_p, y_intc_p, (e.real/np.linalg.norm(e.real))[0][0], (e.real/np.linalg.norm(e.real))[0][1], scale=20, color='b')

                            #                         print('ereal, surfnorm', e.real, line22_normal)
                            e_dot_surf_norm = np.dot(e.real, line22_normal)
                            #                         print("e.surfnorm", e_dot_surf_norm)

                            if e_dot_surf_norm >= 0:
                                #                             print("One satisfies")
                                # calculate new position using 1-dt_frac, u_temp at intersection and x_temp
                                u_emission = line22_normal * np.sqrt(
                                    2 * v_init * q0 / m0)  # <- velocity with which particle is emitted from surface
                                particles.u[ind] = u_emission + q0 / m0 * (
                                            e.real + cross([u_emission], b.real)) * dt * (1 - dt_frac)
                                particles.x[ind] = [x_intc_p, y_intc_p] + particles.u[ind] * dt * (1 - dt_frac)
                            #                             print(particles.x[ind], particles.u[ind])
                            else:
                                lost_particles_indx.append(ind)
            #                             print("|"*50)
            ###################################################################

            #                         surf_pts_neigs = get_neighbours(xsurf, idx)

            #                         # check for intersection of the tentative new point with boundary
            #                         for pt_neigs_ind in range(1, len(surf_pts_neigs)):
            #                             line1 = (particles.x_temp[ind], particles.x_old[ind])
            #                             line2 = (surf_pts_neigs[pt_neigs_ind-1], surf_pts_neigs[pt_neigs_ind])

            #                             bool_intc, (x_intc, y_intc) = segment_intersection(line1, line2)
            #                             if bool_intc:  # <- if true, next position is outside region
            #                                 lost_particles_indx.append(ind)
            #                                 print("|"*50)

            # finally check if particle is at the other boundaries not the wall surface
            for indx_ob, ptx in enumerate(particles.x):
                if ptx[1] < ymin:  # <- bottom edge (rotation axis) check
                    lost_particles_indx.append(indx_ob)
                if ptx[0] < xmin or ptx[0] > xmax:  # <- left and right boundaries
                    lost_particles_indx.append(indx_ob)

            # remove lost points
            if len(lost_particles_indx) != 0:
                particles.remove(lost_particles_indx)

            particles.update_record()
            counter += 1
            t += dt

    print(f'Epk: {epk}, Number of particles left: {particles.len}')
    no_of_remaining_particles.append(particles.len)

# plt.xlim(-.06, .06)
# plt.ylim(0.0, .11)
plt.show()
# print(particles.record)