
import copy
import time
import matplotlib
from matplotlib import cm
import numpy as np


def get_cmap(n, name='jet'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return cm.get_cmap(name, n)


#     return matplotlib.colormaps[name]

def cross(a, b):
    c1 = np.array(a)[:, 1] * np.array(b)[:, 0]
    c2 = -np.array(a)[:, 0] * np.array(b)[:, 0]
    return np.array([c1, c2]).T


def dot(a, b):
    return np.atleast_2d(np.sum(a * b, axis=1)).T


def norm(a):
    return np.atleast_2d(np.linalg.norm(a, axis=1)).T


def segment_intersection(line1, line2):
    x1, y1 = line1[0]
    x2, y2 = line1[1]
    x3, y3 = line2[0][:, 0], line2[0][:, 1]
    x4, y4 = line2[1][:, 0], line2[1][:, 1]

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    #     print('denonm: ', denom)
    #     if denom == 0:
    #        return False, (0, 0)

    tt = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    uu = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / denom

    # get index of where condition is true. condition should be true at only one point
    condition = np.where((tt >= 0) * (tt <= 1) * (uu >= 0) * (uu <= 1))[0]
    #     print('condition', condition)

    if len(condition) == 1:  # <- multiply conditions to combine
        px, py = x1 + tt[condition[0]] * (x2 - x1), y1 + tt[condition[0]] * (y2 - y1)
        return True, np.array([px, py]), condition[0]
    else:
        return False, np.array([0, 0]), 0


def get_neighbours(surf_pts, idx):
    """

    Parameters
    ----------
    surf_pts
    idx

    Returns
    -------

    """
    surf_pts_neigs = np.array(surf_pts[idx])
    #     print("indx", surf_pts_in_box_and_neigs)
    return surf_pts_neigs[surf_pts_neigs[:, 0].argsort()]


def plot_path(particles):
    #     print('before', particles.paths.shape, particles.paths_count)
    particles.paths = np.vstack((particles.paths, particles.x))
    particles.paths_count += 1


#     print('\tafter', particles.paths.shape, particles.paths_count, particles.x.shape)

def hit_bound(particles, t, dt, em, scale):
    #     start = time.time()
    xsurf = particles.bounds

    # check if particle close to boundary
    #     ss = time.time()
    res, indx = particles.distance(100)
    #     print('\t\t Distance: ', time.time() - ss)

    lost_particles_indx = []
    reflected_particles_indx = []
    #     vmag = np.linalg.norm(particles.u, axis=1)
    for ind, (r, idx) in enumerate(zip(res, indx)):
        if r < 5e-2:  # point at boundary, calculate new field value
            # check if point is inside or outside of region
            # get surface points neighbours
            #             ss = time.time()
            surf_pts_neigs = get_neighbours(xsurf, idx)
            #             print(f'\t\t\t {len(surf_pts_neigs)} neigbours: ', time.time() - ss)

            # check for intersection
            # get intersection with old point. loop through points again.
            # the surface edge a line between an outside point and the origin intersects
            # might be different from that with which the line between old and new point intersects

            line11 = (particles.x[ind], particles.x_old[ind])  # <- straight line btw current and previous points
            line22 = surf_pts_neigs[1:], surf_pts_neigs[:-1]
            bool_intc_p, x_intc_p, intc_indx = segment_intersection(line11, line22)

            if bool_intc_p:

                #                 # calculate distance
                #                 kappa = lmbda/(2*np.pi)
                #                 dist = np.sqrt(np.linalg.norm(particles.x_init[ind] - x_intc_p)**2 + kappa*np.linalg.norm(np.exp(1j*particles.phi_init[ind])-np.exp(1j*particles.phi[ind]))**2)
                #                 print('\t\t', dist)
                #                 if True:
                #                    ax.plot([particles.x[ind][0], particles.x_old[ind][0]],
                #                             [particles.x[ind][1], particles.x_old[ind][1]], c='b')
                #                    ax.scatter(x_intc_p[0], x_intc_p[1], color='y', s=100, marker='P')
                #                 calculate time fraction required to hit wall and calculate the field values and
                #                 velocity at this time fraction
                #                 calculate parameter tfrac

                dt_frac = np.linalg.norm(x_intc_p - particles.x_old[ind]) / np.linalg.norm(
                    particles.x[ind] - particles.x_old[ind])
                t_frac = t - dt * (1 - dt_frac)

                #  calculate field values at this time which is a (fraction of dt) + t
                e = scale * np.array([em.e(mesh(*x_intc_p))]) * np.exp(1j * (w * t_frac + particles.phi[ind]))
                b = mu0 * scale * np.array([em.h(mesh(*x_intc_p))]) * np.exp(1j * (w * t_frac + particles.phi[ind]))

                # check if the e-field surface normal is close to zero indicating a possible change in field
                particles.x_temp[ind] = particles.x_old[ind] + particles.u[ind] * dt * dt_frac

                line22 = np.array(line22)[:, intc_indx]
                line22 = line22[line22[:, 0].argsort()]
                line22_normal = -np.array([-(line22[1][1] - line22[0][1]), line22[1][0] - line22[0][0]])

                e_dot_surf_norm = np.dot(e.real, line22_normal)
                if e_dot_surf_norm >= 0:

                    particles.u_temp[ind] = particles.u_old[ind] + q0 / m0 * np.sqrt(
                        1 - (norm([particles.u_old[ind]]) / c0) ** 2) * (
                                                        e.real + cross([particles.u_old[ind]], b.real) - (
                                                            1 / c0 ** 2) * (dot([particles.u_old[ind]], e.real) *
                                                                            particles.u_old[ind])) * dt * dt_frac

                    # check if conditions support secondary electron yield
                    # calculate electron energy
                    umag = np.linalg.norm(particles.u_temp[ind])
                    gamma = 1 / (np.sqrt(1 - (umag / c0) ** 2))
                    pm = gamma * m0 * umag
                    Eq = (gamma - 1) * m0 * c0 ** 2 * 6.241509e18
                    #                     print('\t\t\t\t Energy: ', Eq, " eV", "gamma: ", gamma)
                    #                     if Eq < 50 or Eq > 1500:
                    #     #                     print('\t\t\t\t\tRemoved')
                    #                         lost_particles_indx.append(ind)

                    # calculate new position using 1-dt_frac, u_temp at intersection and x_temp
                    u_emission = line22_normal * np.sqrt(
                        2 * v_init * q0 / m0)  # <- velocity with which particle is emitted from surface
                    #                     particles.u[ind] = u_emission + q0/m0*(e.real + cross([u_emission], b.real))*dt*(1-dt_frac)
                    particles.u[ind] = u_emission + q0 / m0 * np.sqrt(1 - (norm([u_emission]) / c0) ** 2) * (
                                e.real + cross([u_emission], b.real) - (1 / c0 ** 2) * (
                                    dot([u_emission], e.real) * u_emission)) * dt * (1 - dt_frac)
                    particles.x[ind] = x_intc_p + particles.u[ind] * dt * (1 - dt_frac)
                    reflected_particles_indx.append(ind)

                #                     print(f'\t\t\t\t ind: {ind}, nhit: {particles.nhit[ind]+1}')
                else:
                    lost_particles_indx.append(ind)

    # finally check if particle is at the other boundaries not the wall surface
    for indx_ob, ptx in enumerate(particles.x):
        if ptx[1] < ymin:  # <- bottom edge (rotation axis) check
            lost_particles_indx.append(indx_ob)
        if ptx[0] < xmin or ptx[0] > xmax:  # <- left and right boundaries
            lost_particles_indx.append(indx_ob)

    # remove lost points
    #     print(lost_particles_indx, reflected_particles_indx)
    #     print("\thit_bound Exec time: ", time.time() - start)
    return lost_particles_indx, reflected_particles_indx


def collision(active_interval):
    pass


class EMField:
    def __init__(self, e, h):
        #         self.phi = np.atleast_2d(phi).T
        self.e = e
        self.h = h

        # record initial phases


#         self.phi_init = copy.deepcopy(self.phi)

#     def remove(self, ind):
#         self.phi = np.delete(self.phi, ind, axis=0)
#         self.phi_init = np.delete(self.phi_init, ind, axis=0)
#         self.len = len(self.phi)


class Integrators:

    def forward_euler(self, particles, tn, h, em, scale):
        ku1 = h * self.lorentz_force(particles, tn, em, scale)
        particles.u = particles.u + ku1
        particles.x = particles.x + h * particles.u
        plot_path(particles)

        # check for lost particles
        lpi, rpi = hit_bound(particles, tn, h, em, scale)

        if len(rpi) != 0:
            particles.update_hit_count(rpi)

        if len(lpi) != 0:
            particles.remove(lpi)
            em.remove(lpi)

    def implicit_euler(self, particles, tn, h, em, scale):
        pass

    def rk2(self, particles, tn, h, em, scale):

        # k1
        ku1 = h * self.lorentz_force(particles, tn, em, scale)
        kx1 = h * particles.u

        particles_dummy = copy.deepcopy(particles)
        particles_dummy.save_old()
        particles_dummy.u += ku1 / 2
        particles_dummy.x += kx1 / 2

        # check for lost particles
        lpi, rpi = hit_bound(particles_dummy, tn + h / 2, h / 2, em, scale)

        for rp in rpi:
            ku1[rp] = particles_dummy.u[rp] - particles.u[rp]
            kx1[rp] = particles_dummy.x[rp] - particles.x[rp]

        if len(lpi) != 0:
            [ku1], [kx1] = self.rk_update_k([ku1], [kx1], lpi)
            particles_dummy.remove(lpi)
            particles.remove(lpi)
            em.remove(lpi)

        # check if all particles are lost
        if particles.len == 0:
            return False
        # k2=================================================
        ku2 = h * self.lorentz_force(particles_dummy, tn + 2 * h / 3, em,
                                     scale)  # <- particles dummy = particles.u + kn
        kx2 = h * (particles.u + 2 * ku1 / 3)

        particles_dummy = copy.deepcopy(particles)
        particles_dummy.save_old()
        particles_dummy.u += 2 * ku2 / 3
        particles_dummy.x += 2 * kx2 / 3

        # check for lost particles
        lpi, rpi = hit_bound(particles_dummy, tn + 2 * h / 3, 2 * h / 3, em, scale)

        for rp in rpi:
            ku2[rp] = particles_dummy.u[rp] - particles.u[rp]
            kx2[rp] = particles_dummy.x[rp] - particles.x[rp]

        if len(lpi) != 0:
            particles_dummy.remove(lpi)
            [ku1, ku2], [kx1, kx2] = self.rk_update_k([ku1, ku2], [kx1, kx2], lpi)
            particles.remove(lpi)
            em.remove(lpi)

        particles.u = particles.u + ku2
        particles.x = particles.x + kx1

        # check for lost particles
        lpi, rpi = hit_bound(particles, tn, h, em, scale)
        if len(lpi) != 0:
            particles.remove(lpi)
            em.remove(lpi)

        # check if all particles are lost
        if particles.len == 0:
            return False
        plot_path(particles)

    def rk2_23(self, particles, tn, h, em, scale):
        # k1
        ku1 = h * self.lorentz_force(particles, tn, em, scale)
        kx1 = h * particles.u

        particles_dummy = copy.deepcopy(particles)
        particles_dummy.save_old()
        particles_dummy.u += ku1 / 2
        particles_dummy.x += kx1 / 2

        # check for lost particles
        lpi, rpi = hit_bound(particles_dummy, tn + h / 2, h / 2, em, scale)

        for rp in rpi:
            ku1[rp] = particles_dummy.u[rp] - particles.u[rp]
            kx1[rp] = particles_dummy.x[rp] - particles.x[rp]
            particles.update_hit_count(rpi)

        if len(lpi) != 0:
            [ku1], [kx1] = self.rk_update_k([ku1], [kx1], lpi)
            particles_dummy.remove(lpi)
            particles.remove(lpi)
            em.remove(lpi)

        # check if all particles are lost
        if particles.len == 0:
            return False
        # k2=================================================
        ku2 = h * self.lorentz_force(particles_dummy, tn + h / 2, em, scale)  # <- particles dummy = particles.u + kn
        kx2 = h * (particles.u + ku1 / 2)

        particles_dummy = copy.deepcopy(particles)
        particles_dummy.save_old()
        particles_dummy.u += ku2 / 2
        particles_dummy.x += kx2 / 2

        # check for lost particles
        lpi, rpi = hit_bound(particles_dummy, tn + h / 2, h / 2, em, scale)

        for rp in rpi:
            ku2[rp] = particles_dummy.u[rp] - particles.u[rp]
            kx2[rp] = particles_dummy.x[rp] - particles.x[rp]
            particles.update_hit_count(rpi)

        if len(lpi) != 0:
            particles_dummy.remove(lpi)
            [ku1, ku2], [kx1, kx2] = self.rk_update_k([ku1, ku2], [kx1, kx2], lpi)
            particles.remove(lpi)
            em.remove(lpi)

        particles.u = particles.u + (1 / 4 * ku1 + 3 / 4 * ku2)
        particles.x = particles.x + (1 / 4 * kx1 + 3 / 4 * kx2)

        # check for lost particles
        lpi, rpi = hit_bound(particles, tn, h, em, scale)
        if len(lpi) != 0:
            particles.remove(lpi)
            em.remove(lpi)

        # check if all particles are lost
        if particles.len == 0:
            return False
        plot_path(particles)

    def rk4(self, particles, tn, h, em, scale):
        #         start = time.time()
        # k1
        ku1 = h * self.lorentz_force(particles, tn, em, scale)
        kx1 = h * particles.u

        #         ss = time.time()
        particles_dummy = copy.deepcopy(particles)

        particles_dummy.save_old()
        particles_dummy.u += ku1 / 2
        particles_dummy.x += kx1 / 2

        try:
            ku2 = h * self.lorentz_force(particles_dummy, tn + h / 2, em,
                                         scale)  # <- particles dummy = particles.u + kn
            kx2 = h * (particles.u + ku1 / 2)
        except:
            # take full step - euler
            particles.u = particles.u + ku1
            particles.x = particles.x + kx1

            lpi, rpi = hit_bound(particles, tn, h, em, scale)

            if len(rpi) != 0:
                particles.update_hit_count(list(set(rpi)))

            if len(lpi) != 0:
                particles.remove(lpi)

            plot_path(particles)

            return

        particles_dummy = copy.deepcopy(particles)

        particles_dummy.save_old()
        particles_dummy.u += ku2 / 2
        particles_dummy.x += kx2 / 2

        try:
            ku3 = h * self.lorentz_force(particles_dummy, tn + h / 2, em,
                                         scale)  # <- particles dummy = particles.u + kn
            kx3 = h * (particles.u + ku2 / 2)
        except:
            # take full step - euler
            particles.u = particles.u + ku1
            particles.x = particles.x + kx1

            lpi, rpi = hit_bound(particles, tn, h, em, scale)

            if len(rpi) != 0:
                particles.update_hit_count(list(set(rpi)))

            if len(lpi) != 0:
                particles.remove(lpi)

            plot_path(particles)

            return

        particles_dummy = copy.deepcopy(particles)

        particles_dummy.save_old()
        particles_dummy.u += ku3
        particles_dummy.x += kx3

        try:
            ku4 = h * self.lorentz_force(particles_dummy, tn + h, em, scale)  # <- particles dummy = particles.u + kn
            kx4 = h * (particles_dummy.u)
        except:
            # take full step - euler
            particles.u = particles.u + ku1
            particles.x = particles.x + kx1

            lpi, rpi = hit_bound(particles, tn, h, em, scale)

            if len(rpi) != 0:
                particles.update_hit_count(list(set(rpi)))

            if len(lpi) != 0:
                particles.remove(lpi)

            plot_path(particles)

            return

        particles.u = particles.u + 1 / 6 * (ku1 + 2 * ku2 + 2 * ku3 + ku4)
        particles.x = particles.x + 1 / 6 * (kx1 + 2 * kx2 + 2 * kx3 + kx4)

        # check for lost particles
        lpi, rpi = hit_bound(particles, tn, h, em, scale)

        if len(rpi) != 0:
            particles.update_hit_count(list(set(rpi)))

        if len(lpi) != 0:
            particles.remove(lpi)

        # check if all particles are lost
        if particles.len == 0:
            return False
        #         print("rk4 exec time: ", time.time() - start)
        #         print('='*80)

        plot_path(particles)

    def rkf45(self):
        pass

    def rk_update_k(self, ku_list, kx_list, lpi):
        ku_list_new, kx_list_new = [], []
        for ku in ku_list:
            ku = np.delete(ku, lpi, axis=0)
            ku_list_new.append(ku)
        for kx in kx_list:
            kx = np.delete(kx, lpi, axis=0)
            kx_list_new.append(kx)

        return ku_list_new, kx_list_new

    def adams_bashforth(self):
        pass

    def leapfrog(self):
        pass

    def lorentz_force(self, particles, tn, em, scale):
        # get e and b field from eigenmode analysis at particle current position
        ss = time.time()
        #         print(pos, scale)
        e = scale * em.e(mesh(particles.x[:, 0], particles.x[:, 1])) * np.exp(1j * (w * tn + particles.phi))
        b = mu0 * scale * em.h(mesh(particles.x[:, 0], particles.x[:, 1])) * np.exp(1j * (w * tn + particles.phi))

        #         k = q0/m0*np.sqrt(1 - (norm(u)/c0)**2)*(e.real + cross(u, b.real)-(1/c0**2)*(dot(u, e.real)*u))  # <- relativistic
        k = q0 / m0 * np.sqrt(1 - (norm(particles.u) / c0) ** 2) * (
                    e.real + cross(particles.u, b.real) - (1 / (c0 ** 2)) * (
                        dot(particles.u, e.real) * particles.u))  # <- relativistic
        #         print('\t\t lorentz force: ', time.time() - ss, len(pos))
        return k


class Particles:
    def __init__(self, xrange, init_v, bounds, phi, cmap='jet'):
        self.cmap = cmap
        M = len(phi)

        self.bounds = np.array(bounds)

        self.x = self.bounds[(self.bounds[:, 0] > xrange[0]) & (self.bounds[:, 0] < xrange[1])]
        shape = self.x.shape
        self.len = len(self.x)

        # get normal pointing inwards of emission points
        self.pt_normals = np.ones(self.x.shape)
        # get neighbouring points
        res, idxs = self.distance(1)
        for nn, idx in enumerate(idxs):
            n12 = self.get_point_normal(idx)
            self.pt_normals[nn] = n12 / np.linalg.norm(n12)

        # repeat into multidomensional array
        self.x = np.array(self.x.tolist() * M)
        self.pt_normals = np.array(self.pt_normals.tolist() * M)

        # convert velocity from eV to m/s
        init_v = np.sqrt(2 * init_v * q0 / m0)
        self.u = self.pt_normals * init_v

        self.phi = np.atleast_2d(np.repeat(phi, self.len)).T

        cmap = get_cmap(len(self.x), self.cmap)
        self.colors = np.array([cmap(i) for i in range(len(self.x))])

        #         self.len = self.x[M,S,:]

        self.x_old, self.u_old, self.phi_old = np.zeros(self.x.shape), np.zeros(self.u.shape), np.zeros(
            self.phi.shape)  # <- hold solution from prev time step
        self.x_temp, self.u_temp, self.phi_temp = np.zeros(self.x.shape), np.zeros(self.u.shape), np.zeros(
            self.phi.shape)  # <- hold tentative change in position

        # record initial particles state
        self.x_init = copy.deepcopy(self.x)
        self.u_init = copy.deepcopy(self.u)
        self.phi_init = copy.deepcopy(self.phi)

        # particles path
        self.paths = copy.deepcopy(self.x)
        self.paths_count = 1

        self.record = [self.x]
        self.lost_particles = []
        self.nhit = np.zeros(len(self.x))

        self.bright_set = []
        self.shadow_set = []

    def save_old(self):
        self.x_old = copy.deepcopy(self.x)
        self.u_old = copy.deepcopy(self.u)
        self.phi_old = copy.deepcopy(self.phi)

        self.x_temp = copy.deepcopy(self.x)
        self.u_temp = copy.deepcopy(self.u)
        self.phi_temp = copy.deepcopy(self.phi)

    def distance(self, n):
        ss = time.time()
        #         closest_dist_to_surface1 = []
        #         indx1 = []
        #         for xx in self.x:
        #             res = np.linalg.norm(xx - self.bounds, axis=1)
        #             closest_dist_to_surface1.append(np.min(res))
        #             indx1.append([list(np.argpartition(res, n)[0:n])][0])  # <- index of two closest surface points
        #         print(time.time()-ss)
        # #         print(np.sort(resl))
        #         print(np.atleast_2d(closest_dist_to_surface).T)
        #         print(type(indx1))
        #         print()
        #         ss = time.time()

        norms = np.linalg.norm(self.x[:, None] - self.bounds, axis=-1)
        closest_dist_to_surface = norms.min(axis=1)
        indx = norms.argsort(axis=1).T[0:n, :]
        #         print('v', time.time()-ss)
        #         print()
        #         print('\t vectorised', np.sort(norms))
        #         print('\t vectorised', np.atleast_2d(closest_dist_to_surface).T)
        #         print('\t vectorised')
        #         print(indx.T.tolist())
        #         print(type(indx))
        #         print()
        return np.atleast_2d(closest_dist_to_surface).T, indx.T.tolist()

    def distance2(self, n):
        ss = time.time()
        closest_dist_to_surface1 = []
        indx1 = []
        for xx, xx_old in zip(self.x, self.x_old):
            box_bound = self.bounds[(self.bounds[:, 0] > min(xx[0], xx_old[0]) - 0.005) &
                                    (self.bounds[:, 0] < max(xx[0], xx_old[0]) + 0.005) &
                                    (self.bounds[:, 1] > min(xx[1], xx_old[1]) - 0.01) &
                                    (self.bounds[:, 1] < max(xx[1], xx_old[1]) + 0.01)]
            #             print(len(box_bound), xx, xx_old)
            res = np.linalg.norm(xx - box_bound, axis=1)

        #             closest_dist_to_surface1.append(np.min(res))
        #             indx1.append([list(np.argpartition(res, n)[0:n])][0])  # <- index of two closest surface points
        print("bound", time.time() - ss)
        print()
        print()

        return np.atleast_2d(closest_dist_to_surface).T, indx.T.tolist()

    def remove(self, ind, bright='no'):
        if bright != 'yes':
            # add to shadow set before removal from main set
            self.shadow_set.append(self.paths[[ii * len(self.x) + np.array(ind) for ii in range(self.paths_count)]])

        self.paths = np.delete(self.paths, [ii * len(self.x) + np.array(ind) for ii in range(self.paths_count)], axis=0)

        self.x = np.delete(self.x, ind, axis=0)
        self.u = np.delete(self.u, ind, axis=0)

        self.x_old = np.delete(self.x_old, ind, axis=0)
        self.u_old = np.delete(self.u_old, ind, axis=0)

        self.phi = np.delete(self.phi, ind, axis=0)

        self.x_temp = np.delete(self.x_temp, ind, axis=0)
        self.u_temp = np.delete(self.u_temp, ind, axis=0)

        self.x_init = np.delete(self.x_init, ind, axis=0)
        self.u_init = np.delete(self.u_init, ind, axis=0)
        self.phi_init = np.delete(self.phi_init, ind, axis=0)

        self.colors = np.delete(self.colors, ind, axis=0)

        self.len = len(self.x)

        # print number of hits of particle before deleting
        # print('\t\t\t\t\t', self.nhit[ind])
        self.nhit = np.delete(self.nhit, ind, axis=0)

    def colors(self):
        return self.colors

    def set_cmap(self, cmap):
        self.cmap = cmap
        cmap = get_cmap(len(self.x), cmap)
        self.colors = np.array([cmap(i) for i in range(len(self.x))])

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

    def update_hit_count(self, inds):
        self.nhit[inds] = self.nhit[inds] + 1

        # sort inds to start deleting from largest index
        inds.sort(reverse=True)
        # check if nhit = 20 and remove from main set and add to bright set
        for ind in inds:
            if self.nhit[ind] == 20:
                self.bright_set.append(self.paths[[ii * len(self.x) + np.array(ind) for ii in range(self.paths_count)]])
                # remove the index from main set
                self.remove(ind, bright='yes')


n = 1
q0 = 1.60217663e-19
m0 = 9.1093837e-31
mu0 = 4 * pi * 1e-7
eps0 = 8.85418782e-12
c0 = 299792458

# get surface points
pec_boundary = mesh.Boundaries("default")
bel = [xx.vertices for xx in pec_boundary.Elements()]
bel_unique = list(set(itertools.chain(*bel)))
xpnts_surf = sorted([mesh.vertices[xy.nr].point for xy in bel_unique])
xsurf = np.array(xpnts_surf)

# get xmin, xmax, ymin
xmin = face.vertices.Min(X).p[0]
xmax = face.vertices.Max(X).p[0]
ymin = 0

# define one internal point
xin = [0, 0]

integrators = Integrators()

no_of_remaining_particles = []
# epks_v = 1/Epk*1e6*np.linspace(1, 60, 119)
# phi_v = np.linspace(0, 2*np.pi, 72)  # <- initial phase
epks_v = 1 / Epk * 1e6 * np.linspace(30, 60, 1)
phi_v = np.linspace(0, 2 * np.pi, 72)  # <- initial phase

v_init = 2  # eV

# calculate time for 10 cycles, 20 alternations
T = 1 / (freq_fes[n] * 1e6) * 10
lmbda = c0 / (freq_fes[n] * 1e6)
w = 2 * np.pi * freq_fes[n] * 1e6

particles_left = []
particles_nhits = []
for epk in epks_v:
    scale = epk  # <- scale Epk to 1 MV/m and multiply by sweep value
    t = 0
    dt = 1e-11
    PLOT = False

    e_check_list = []
    time_list = []
    error = False
    counter = 0

    init_pos = [-0.005, 0.005]
    particles = Particles(init_pos, v_init, xsurf, phi_v, cmap='jet')

    print('Number of particles: ', len(particles.x))

    em = EMField(copy.deepcopy(gfu_E[n]), copy.deepcopy(gfu_H[n]))

    # move particles with initial velocity. ensure all initial positions after first move lie inside the bounds
    particles.x = particles.x + particles.u * dt
    # particles2.x = particles2.x + particles2.u*dt

    record = {}

    #     while t < T:
    while t < 1000e-10:  # particles.len > 0:
        #         start = time.time()
        if particles.len != 0:
            particles.save_old()
            #             ss = time.time()
            integrators.rk4(particles, t, dt, em, scale)
            #             print('\t\t\trk: ', time.time()-ss)
            particles.update_record()

        #         if particles.len != 0:
        #             particles.save_old()
        #             integrators.forward_euler(particles, t, dt, em, scale)
        #             particles.update_record()

        # update record

        #         # define dataframe to hold values
        #         dd = {"scale": scale,
        #               "phi": particles.phi[:, 0],
        #               "p.x0": particles.x[:, 0],
        #               "p.x1": particles.x[:, 1],
        #               "p.u0": particles.u[:, 0],
        #               "p.u1": particles.u[:, 1]}

        #         df = pd.DataFrame(dd)
        #         record[t] = df

        #     print('One loop time: ', time.time()-start)
        #         print('nhits: ', particles.nhit)
        counter += 1
        t += dt

    if len(particles.nhit) == 0:
        particles_nhits.append(0)
    else:
        particles_nhits.append(particles.nhit[0])

    particles_left.append(len(particles.x))
    print(f"Epk: {epk}, particles left: {len(particles.x)}")

# print(record)
# print(df[(df['phi'] == 0) & (df['scale'] == em.scale[0, 0])])

print("Done")