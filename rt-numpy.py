import numpy as np

from math_tools import normalize
from rng import RNG
from sampling import cosine_weighted_sample_on_hemisphere
from sphere import Sphere
from specular import ideal_specular_reflect, ideal_specular_transmit

# image I/O
def write_ppm(w, h, Ls, fname = "numpy-image.ppm"):
    with open(fname, 'w') as outfile:
        outfile.write('P3\n{0} {1}\n{2}\n'.format(w, h, 255));
        for i in range(Ls.shape[0]):
            outfile.write('{0} {1} {2} '.format(to_byte(Ls[i,0]), to_byte(Ls[i,1]), to_byte(Ls[i,2])));

# ray bundles
class Rays(object):

    def __init__(self, Os, Ds):
        self.Os = np.copy(Os)
        self.Ds = np.copy(Ds)

    def __call__(self, t):
        return self.Os + self.Ds * t

    def __str__(self):
        return 'o: ' + str(self.o) + '\n' + 'd: ' + str(self.d) + '\n'



class Sphere(object):

    EPSILON_SPHERE = 1e-4

    def __init__(self, r, p, e = np.zeros((3), dtype=np.float64), f = np.zeros((3), dtype=np.float64)):
        self.r = np.float64(r)
        self.p = np.copy(p)
        self.e = np.copy(e)
        self.f = np.copy(f)

    def intersect(self, rays):
        # (o + t*d - p) . (o + t*d - p) - r*r = 0
        # <=> (d . d) * t^2 + 2 * d . (o - p) * t + (o - p) . (o - p) - r*r = 0
        # 
        # Discriminant check
        # (2 * d . (o - p))^2 - 4 * (d . d) * ((o - p) . (o - p) - r*r) <? 0
        # <=> (d . (o - p))^2 - (d . d) * ((o - p) . (o - p) - r*r) <? 0
        # <=> (d . op)^2 - 1 * (op . op - r*r) <? 0
        # <=> b^2 - (op . op) + r*r <? 0
        # <=> D <? 0
        #
        # Solutions
        # t = (- 2 * d . (o - p) +- 2 * sqrt(D)) / (2 * (d . d))
        # <=> t = dop +- sqrt(D)

        op = self.p - rays.Os
        dop = rays.Ds.dot(op)
        D = dop * dop - op.dot(op) + self.r * self.r

        breakpoint()

        return False

        if D < 0:
            return False

        sqrtD = np.sqrt(D)

        tmin = dop - sqrtD
        if (ray.tmin < tmin and tmin < ray.tmax):
            ray.tmax = tmin
            return True

        tmax = dop + sqrtD
        if (ray.tmin < tmax and tmax < ray.tmax):
            ray.tmax = tmax
            return True
        
        return False



# Scene
REFRACTIVE_INDEX_OUT = 1.0
REFRACTIVE_INDEX_IN = 1.5

spheres = [
        #Sphere(1e5,  np.array([1e5 + 1, 40.8, 81.6],    dtype=np.float64), f=np.array([0.75,0.25,0.25],      dtype=np.float64)),
	    #Sphere(1e5,  np.array([-1e5 + 99, 40.8, 81.6],  dtype=np.float64), f=np.array([0.25,0.25,0.75],      dtype=np.float64)),
	    #Sphere(1e5,  np.array([50, 40.8, 1e5],          dtype=np.float64), f=np.array([0.75, 0.75, 0.75],    dtype=np.float64)),
	    #Sphere(1e5,  np.array([50, 40.8, -1e5 + 170],   dtype=np.float64)),
	    #Sphere(1e5,  np.array([50, 1e5, 81.6],          dtype=np.float64), f=np.array([0.75, 0.75, 0.75],    dtype=np.float64)),
	    #Sphere(1e5,  np.array([50, -1e5 + 81.6, 81.6],  dtype=np.float64), f=np.array([0.75, 0.75, 0.75],    dtype=np.float64)),
	    #Sphere(16.5, np.array([27, 16.5, 47],           dtype=np.float64), f=np.array([0.999, 0.999, 0.999], dtype=np.float64), reflection_t=Sphere.Reflection_t.SPECULAR),
	    #Sphere(16.5, np.array([73, 16.5, 78],           dtype=np.float64), f=np.array([0.999, 0.999, 0.999], dtype=np.float64), reflection_t=Sphere.Reflection_t.REFRACTIVE),
	    Sphere(600,  np.array([50, 681.6 - .27, 81.6],  dtype=np.float64), e=np.array([12, 12, 12],          dtype=np.float64))
        ]


def intersect(ray):
    id = None
    hit = False
    for i in range(len(spheres)):
        if spheres[i].intersect(ray):
            hit = True
            id = i
    return hit, id


def radiance(ray, rng):
    r = ray
    L = np.zeros((3), dtype=np.float64)
    F = np.ones((3), dtype=np.float64)

    hit, id = intersect(r)
    if(not hit):
        return L

    return np.array([1,0,0])
    
    if False: # todo
        while (True):
            hit, id = intersect(r)
            if (not hit):
                return L

            shape = spheres[id]
            p = r(r.tmax)
            n = normalize(p - shape.p)

            L += F * shape.e
            F *= shape.f
            
    	    # Russian roulette
            if r.depth > 4:
                continue_probability = np.amax(shape.f)
                if rng.uniform_float() >= continue_probability:
                    return L
                F /= continue_probability

            # Next path segment
            if shape.reflection_t == Sphere.Reflection_t.SPECULAR:
                d = ideal_specular_reflect(r.d, n)
                r = Ray(p, d, tmin=Sphere.EPSILON_SPHERE, depth=r.depth + 1)
                continue
            elif shape.reflection_t == Sphere.Reflection_t.REFRACTIVE:
                d, pr = ideal_specular_transmit(r.d, n, REFRACTIVE_INDEX_OUT, REFRACTIVE_INDEX_IN, rng)
                F *= pr
                r = Ray(p, d, tmin=Sphere.EPSILON_SPHERE, depth=r.depth + 1)
                continue
            else:
                w = n if n.dot(r.d) < 0 else -n
                u = normalize(np.cross(np.array([0.0, 1.0, 0.0], np.float64) if np.fabs(w[0]) > 0.1 else np.array([1.0, 0.0, 0.0], np.float64), w))
                v = np.cross(w, u)

                sample_d = cosine_weighted_sample_on_hemisphere(rng.uniform_float(), rng.uniform_float())
                d = normalize(sample_d[0] * u + sample_d[1] * v + sample_d[2] * w)
                r = Ray(p, d, tmin=Sphere.EPSILON_SPHERE, depth=r.depth + 1)
                continue

import sys

if __name__ == "__main__":
    rng = RNG()
    nb_samples = 1 #int(sys.argv[1]) // 4 if len(sys.argv) > 1 else 1

    w = 5 #1024
    h = 3 #768

    eye = np.array([50, 52, 295.6], dtype=np.float64)
    gaze = normalize(np.array([0, -0.042612, -1], dtype=np.float64))
    fov = 0.5135
    cx = np.array([w * fov / h, 0.0, 0.0], dtype=np.float64)
    cy = normalize(np.cross(cx, gaze)) * fov

    Ls = np.zeros((w * h, 3), dtype=np.float64)
    

    pixel_origins = np.zeros((w * h, 3))
    pixel_directions = np.zeros((w * h, 3))
    

    for y in range(h):
        # pixel row
        print('\rRendering ({0} spp) {1:0.2f}%'.format(nb_samples, 100.0 * y / (h - 1)))
        for x in range(w):
            # pixel column
            #i = (h - 1 - y) * w + x
            i = y * w + x
            L = np.zeros((3), dtype=np.float64)
            for s in range(nb_samples):
                #  samples per subpixel
                u1 = 1.0 #2.0 * rng.uniform_float()
                u2 = 1.0 #2.0 * rng.uniform_float()
                dx = 0.0 #np.sqrt(u1) - 1.0 if u1 < 1 else 1.0 - np.sqrt(2.0 - u1)
                dy = 0.0 #np.sqrt(u2) - 1.0 if u2 < 1 else 1.0 - np.sqrt(2.0 - u2)
                d = cx * (((0 + 0.5 + dx) / 2.0 + x) / w - 0.5) + \
                    cy * (((0 + 0.5 + dy) / 2.0 + y) / h - 0.5) + gaze
                    
                print(i)
                pixel_origins[i, :] = eye + d * 130
                pixel_directions[i, :] = normalize(d)
                #L += radiance(Ray(eye + d * 130, normalize(d), tmin=Sphere.EPSILON_SPHERE), rng) * (1.0 / nb_samples)
            #Ls[i,:] += 0.25 * np.clip(L, a_min=0.0, a_max=1.0)

    eye_rays = Rays(pixel_origins, pixel_directions)
    
    i = np.arange( w * h )
    directions = cx * (((0 + 0.5 + 0.0) / 2.0 + (i[:,np.newaxis]%w)) / w - 0.5) + \
                 cy * (((0 + 0.5 + 0.0) / 2.0 + np.floor(i[:,np.newaxis] / w) ) / h - 0.5) + gaze
    
    vectorized_origins = eye + directions * 130
    #breakpoint()
    
    norms = np.linalg.norm(directions, axis=1)
    normed_directions = directions / norms[:,np.newaxis]
    vectorized_directions = np.where( np.isclose( norms, 0)[:,np.newaxis], directions, normed_directions )
    vectorized_eye_rays = Rays(vectorized_origins, vectorized_directions)
    breakpoint()

    write_ppm(w, h, Ls)
