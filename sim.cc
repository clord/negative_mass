#include <Cocoa/Cocoa.h>
#include <dispatch/dispatch.h>
#import <QTKit/QTKit.h>
#include <vector>
#include <iostream>

const double C = 1.0;
const int num_masses = 1800;


// A force is a vector in 3-space.
struct vec_t {
	double x, y, z;
	vec_t(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
	bool operator!=(const vec_t & v) const {
		return x != v.x ||
		       y != v.y ||
		       z != v.z;
	}
	bool operator==(const vec_t & v) const {
		return x == v.x &&
		       y == v.y &&
		       z == v.z;
	}
	vec_t operator+(const vec_t & v) {
		vec_t r(*this);
		r.x += v.x;
		r.y += v.y;
		r.z += v.z;
		return r;
	}
	vec_t& operator*=(const vec_t & v) {
		x *= v.x;
		y *= v.y;
		z *= v.z;
		return *this;
	}
	vec_t& operator+=(const vec_t & v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	vec_t& operator-=(const vec_t & v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}
	vec_t operator+(const double & s) {
		vec_t r(*this);
		r.x += s;
		r.y += s;
		r.z += s;
		return r;
	}
	vec_t operator*(const double & s) {
		vec_t r(*this);
		r.x *= s;
		r.y *= s;
		r.z *= s;
		return r;
	}
	vec_t operator*(const vec_t & v) {
		vec_t r(*this);
		r.x *= v.x;
		r.y *= v.y;
		r.z *= v.z;
		return r;
	}
	vec_t operator/(const vec_t & v) {
		vec_t r(*this);
		r.x /= v.x;
		r.y /= v.y;
		r.z /= v.z;
		return r;
	}
	vec_t operator/(const double & s) {
		vec_t r(*this);
		r.x /= s;
		r.y /= s;
		r.z /= s;
		return r;
	}
	vec_t operator-(const vec_t & v) {
		vec_t r(*this);
		r.x -= v.x;
		r.y -= v.y;
		r.z -= v.z;
		return r;
	}

	vec_t unit_vector() const {
		vec_t unit(*this);
		double den(sqrt(unit.length_squared()));
		if (den > 0.0000000001) {
			unit.x /= den;
			unit.y /= den;
			unit.z /= den;
		}
		return unit;
	}
	
	double length_squared() const {
		vec_t vec(*this);
		vec = vec * vec;
		return vec.x + vec.y + vec.z;
	}
	
	double distance_squared_to(const vec_t & to) const {
		vec_t minus(*this);
		minus -= to;
		vec_t prod(minus * minus);
		return prod.x + prod.y + prod.z;
	}
};

struct mass_t {
	double mass;
	vec_t position;
	vec_t velocity;
	mass_t() : position(0,0,0), mass(0), velocity(0,0,0) {}
	mass_t(double _x, double _y, double _z, double _t, double _u, double _v, double _m) : 
	                  position(_x, _y, _z), mass(_m), velocity(_t, _u, _v) {}
	bool operator!=(const mass_t & v) const {
		return mass != v.mass || position != v.position || velocity != v.velocity;
	}
	bool operator==(const mass_t & v) const {
		return mass == v.mass && position == v.position && velocity == v.velocity;
	}
};


// Assumes Δt is 1, eliminating some multiplication. 
mass_t translate(const std::vector<mass_t> & masses, const int particle_i)
{
	mass_t particle(masses[particle_i]);

	vec_t vel_vec(0.0, 0.0, 0.0);
	vec_t pos_vec(0.0, 0.0, 0.0);
	for (std::vector<mass_t>::const_iterator i = masses.begin(); i != masses.end(); ++i) {
		mass_t mass_i(*i);
		vec_t vel_diff = particle.velocity - mass_i.velocity;
		vec_t pos_diff = particle.position - mass_i.position;
		
		double dsq = pos_diff.length_squared();
		if (dsq < 0.00000001) continue;
		if (dsq > DBL_MAX) continue;
		double same_dir = ((particle.mass < 0 && mass_i.mass < 0) ||
		                   (particle.mass > 0 && mass_i.mass > 0)) ? -1.0 : 1.0;
		
		double force(((fabs(particle.mass) * fabs(mass_i.mass)) / (dsq)));
		double accel(same_dir * force / fabs(particle.mass));
		
		pos_vec = pos_diff.unit_vector() * (pos_vec + vel_vec + 0.5 * accel);
		vel_vec = vel_diff.unit_vector() * (vel_vec + accel);

	}

	particle.velocity = particle.velocity + vel_vec;
	particle.position = particle.position + pos_vec;
	return particle;
}



std::vector<mass_t> b_ms;
std::vector<mass_t> a_ms;

std::vector<mass_t> * prev_masses = &a_ms;
std::vector<mass_t> * next_masses = &b_ms;
std::vector<mass_t> *tmp;

int main(int argc, char ** argv) {

	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	QTMovie * movie = [[QTMovie alloc] initToWritableFile: @"movie.mov" error: nil];
	[movie retain];
	dispatch_queue_t queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
	srand(time(0));
	// Fill the world with masses
	for (int i = 0; i < num_masses; ++i) {
		double mss = ((rand()/(double)RAND_MAX) - 0.5) * 500.0;
		mass_t m(145.0 + 324.0 * (rand()/(double)RAND_MAX), 
		         90.0 + 324.0 * (rand()/(double)RAND_MAX), 
		         145.0 + 324.0 * (rand()/(double)RAND_MAX),
		         (rand()/(double)RAND_MAX)/20.0, 
		         (rand()/(double)RAND_MAX)/20.0, 
		         (rand()/(double)RAND_MAX)/20.0, 
		         mss);

		prev_masses->push_back(m);
		next_masses->push_back(m);
	}

	QTTime frameDuration = QTMakeTime(1, 60);
	NSDictionary* movieAttributes = [NSDictionary dictionaryWithObjectsAndKeys: 
		@"tiff", QTAddImageCodecType, 
		[NSNumber numberWithLong: codecHighQuality], QTAddImageCodecQuality, nil];

	// Walk through the array, flipping next/prev at each step
	for (int iter = 1; iter <= 23000; ++iter) {
		dispatch_group_t iteration_group(dispatch_group_create());

		std::vector<mass_t> & nm = *next_masses;
		std::vector<mass_t> & pm = *prev_masses;

		const int skip = 20;
		if (iter % skip == 0) {
			dispatch_group_async(iteration_group, queue, ^{
				std::cout << "writing frame " << iter / skip << std::endl;
				NSImage * frame = [[NSImage alloc] initWithSize: NSMakeSize(640.0, 480.0)];
				[frame lockFocus];
			
				[[NSColor colorWithDeviceWhite:0.1 alpha:1.0] set];
				NSRectFill(NSMakeRect(0.0, 0.0, 640.0, 480.0));
				for (std::vector<mass_t>::const_iterator i = pm.begin(); i !=pm.end(); ++i) {
					double scale = i->position.z / 1024.0;
					if (i->mass < 0)
						[[NSColor blackColor] set];
					else
						[[NSColor colorWithDeviceHue: (0.5+i->mass/4000.0) saturation:0.4
						                  brightness:scale+fmax(rand()/(double)RAND_MAX,0.5)/1.6 alpha:1.0] set];
					double rad = 3.5 * scale;
					NSPoint p = {i->position.x - rad, i->position.y - rad};
					NSRect r = {p.x, p.y, rad*2, rad*2};
					NSBezierPath * pat = [NSBezierPath bezierPathWithOvalInRect:r];
					[pat fill];
				}
			
				[frame unlockFocus];
				[movie addImage: frame forDuration: frameDuration withAttributes: movieAttributes];
				[movie updateMovieFile];
				[frame release];
			});
		}
		
		const int work_units = 20;
		for (int s = 0; s < work_units; ++s) {
			dispatch_group_async(iteration_group, queue, ^{
				for (int i = s * num_masses/work_units; i < (s + 1) * num_masses/work_units; ++i) {
					nm[i] = translate(pm, i);
				}
			});
		}
		
		dispatch_group_wait(iteration_group, DISPATCH_TIME_FOREVER);
		dispatch_release(iteration_group);
		//std::cout << "iteration " << iter << " done." << std::endl;
		tmp = prev_masses;
		prev_masses = next_masses;
		next_masses = tmp;
	}
	[movie release];
	[pool release];
	std::cout << "done" << std::endl;
}
