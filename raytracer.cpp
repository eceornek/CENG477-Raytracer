#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <math.h>
#include <limits>

typedef unsigned char RGB[3];

using namespace parser;
using namespace std;

float findLength(const Vec3f &a)
{
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

float findDistance(const Vec3f &a, const Vec3f &b)
{
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}
float dot(const Vec3f &a, const Vec3f &b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3f multS(const Vec3f &a, float scalar)
{
	Vec3f result;
	result.x = a.x * scalar;
	result.y = a.y * scalar;
	result.z = a.z * scalar;
	return result;
}

Vec3f normalizeVector(const Vec3f &a)
{
	return multS(a, (1.0 / sqrt(dot(a, a))));
}

Vec3f subtract(const Vec3f &a, const Vec3f &b)
{
	Vec3f result;
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;
	return result;
}

Vec3f add(const Vec3f &a, const Vec3f &b)
{
	Vec3f result;
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;
	return result;
}

Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
	Vec3f result;
	result.x = a.y * b.z - a.z * b.y;
	result.y = a.z * b.x - a.x * b.z;
	result.z = a.x * b.y - a.y * b.x;
	return result;
}

float determinant(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2)
{
	return (v0.x * (v1.y * v2.z - v2.y * v1.z)) + (v0.y * (v2.x * v1.z - v1.x * v2.z)) + (v0.z * (v1.x * v2.y - v1.y * v2.x));
}

class Ray
{
	Vec3f e;
	Vec3f d;
	Vec3f gaze;
	float minT = numeric_limits<float>::max();
	bool sphere = false;
	bool triangle = false;
	bool mesh = false;
	int objID = -1;
	Vec3f intersectionPoint;
	Vec3f surfaceNormal;
	int material_id;

public:
	Ray(const Camera &camera, int i, int j)
	{
		float left = camera.near_plane.x;
		float right = camera.near_plane.y;
		float bottom = camera.near_plane.z;
		float top = camera.near_plane.w;

		this->gaze = normalizeVector(camera.gaze);

		float su = (right - left) * (j + 0.5) / camera.image_width;
		float sv = (top - bottom) * (i + 0.5) / camera.image_height;

		Vec3f m, q, u, v;
		m.x = camera.position.x + (gaze.x * camera.near_distance);
		m.y = camera.position.y + (gaze.y * camera.near_distance);
		m.z = camera.position.z + (gaze.z * camera.near_distance);

		u = crossProduct(gaze, camera.up);
		u = normalizeVector(u);

		v = crossProduct(u, gaze);

		q.x = m.x + (u.x * left) + (v.x * top);
		q.y = m.y + (u.y * left) + (v.y * top);
		q.z = m.z + (u.z * left) + (v.z * top);

		Vec3f s;
		s.x = q.x + (u.x * su) - (v.x * sv);
		s.y = q.y + (u.y * su) - (v.y * sv);
		s.z = q.z + (u.z * su) - (v.z * sv);

		this->e = camera.position;
		this->d = subtract(s, camera.position);
		this->d = normalizeVector(this->d);
	}
	Ray(Vec3f e, Vec3f d)
	{
		this->e = e;
		this->d = d;
	}

	Vec3f findIntersectionPoint(float t)
	{
		Vec3f result;
		result.x = this->e.x + (t * this->d.x);
		result.y = this->e.y + (t * this->d.y);
		result.z = this->e.z + (t * this->d.z);
		return result;
	}
	void sphereIntersection(const Sphere &sphere, int sphereIndex, Vec3f &center)
	{
		int material_id = sphere.material_id;
		float radius = sphere.radius;

		float A = dot(this->d, this->d);
		Vec3f e_minus_c = subtract(this->e, center);
		float B = dot(this->d, e_minus_c) * 2;
		float C = dot(e_minus_c, e_minus_c) - (radius * radius);
		float discriminant = (B * B) - (4 * A * C);

		if (discriminant >= 0)
		{
			float sqrt_discriminant = sqrtf(discriminant);
			float t1 = ((-1 * B) + sqrt_discriminant) / (2 * A);
			float t2 = ((-1 * B) - sqrt_discriminant) / (2 * A);

			float t = fmin(t1, t2);
			if (t < this->minT && t >= 0)
			{
				this->intersectionPoint = findIntersectionPoint(t);
				this->surfaceNormal = normalizeVector(subtract(this->intersectionPoint, center));
				// this->surfaceNormal.x /= radius;
				// this->surfaceNormal.y /= radius;
				// this->surfaceNormal.z /= radius;
				this->minT = t;
				this->triangle = false;
				this->sphere = true;
				this->mesh = false;
				this->objID = sphereIndex;
				this->material_id = material_id;
			}
		}
	}
	void triangleIntersection(const Triangle &triangle, const Vec3f &a, const Vec3f &b, const Vec3f &c, int triangleIndex)
	{
		int material_id = triangle.material_id;

		Vec3f A_1 = subtract(a, b);
		Vec3f A_2 = subtract(a, c);
		Vec3f A_3 = this->d;

		Vec3f a_minus_e = subtract(a, this->e);

		float detA = determinant(A_1, A_2, A_3);
		if (detA == 0.0)
		{
			return;
		}

		float gamma = (determinant(A_1, a_minus_e, A_3)) / detA;
		float beta = (determinant(a_minus_e, A_2, A_3)) / detA;
		float t = (determinant(A_1, A_2, a_minus_e)) / detA;

		if (gamma < 0 || beta < 0)
		{
			return;
		}
		if (beta > (1 - gamma) || gamma > (1 - beta))
		{
			return;
		}
		if (t > 0 && t < this->minT)
		{
			this->triangle = true;
			this->sphere = false;
			this->mesh = false;
			this->minT = t;
			this->objID = triangleIndex;
			this->material_id = material_id;
			this->intersectionPoint = findIntersectionPoint(t);
			this->surfaceNormal = normalizeVector(crossProduct(subtract(b, a), subtract(c, a)));
		}
	}

	void triangleIntersectionForMesh(const Vec3f &a, const Vec3f &b, const Vec3f &c, int material_id, int meshIndex)
	{

		Vec3f A_1 = subtract(a, b);
		Vec3f A_2 = subtract(a, c);
		Vec3f A_3 = this->d;
		Vec3f a_minus_e = subtract(a, this->e);

		float detA = determinant(A_1, A_2, A_3);
		if (detA == 0.0)
		{
			return;
		}

		float gamma = (determinant(A_1, a_minus_e, A_3)) / detA;
		float beta = (determinant(a_minus_e, A_2, A_3)) / detA;
		float t = (determinant(A_1, A_2, a_minus_e)) / detA;

		if (gamma < 0 || beta < 0)
		{
			return;
		}
		if (beta > (1 - gamma) || gamma > (1 - beta))
		{
			return;
		}

		if (t > 0 && t < this->minT)
		{
			this->triangle = false;
			this->sphere = false;
			this->mesh = true;
			this->minT = t;
			this->objID = meshIndex;
			this->material_id = material_id;
			this->intersectionPoint = findIntersectionPoint(t);
			this->surfaceNormal = normalizeVector(crossProduct(subtract(b, a), subtract(c, a)));
		}
	}

	void meshIntersection(const Mesh &mesh, const Scene &scene, int meshIndex)
	{
		int size = mesh.faces.size();
		for (int faceNumber = 0; faceNumber < size; faceNumber++)
		{
			Vec3f v0 = scene.vertex_data[mesh.faces[faceNumber].v0_id - 1];
			Vec3f v1 = scene.vertex_data[mesh.faces[faceNumber].v1_id - 1];
			Vec3f v2 = scene.vertex_data[mesh.faces[faceNumber].v2_id - 1];

			triangleIntersectionForMesh(v0, v1, v2, mesh.material_id, meshIndex);
		}
	}
	Vec3f computeAmbient(const Scene &scene)
	{
		Vec3f ambient;
		ambient.x = scene.materials[this->material_id - 1].ambient.x * scene.ambient_light.x;
		ambient.y = scene.materials[this->material_id - 1].ambient.y * scene.ambient_light.y;
		ambient.z = scene.materials[this->material_id - 1].ambient.z * scene.ambient_light.z;
		return ambient;
	}

	Vec3f computeColor(const Scene &scene, const Camera &currentCamera, int maxDepth)
	{
		int numberOfLights = scene.point_lights.size();
		int numberOfSpheres = scene.spheres.size();
		int numberOfTriangles = scene.triangles.size();
		int numberOfMeshes = scene.meshes.size();

		float r;
		float g;
		float b;

		Vec3f color;

		if (this->objID != -1) // hit happened
		{
			// Compute ambient only once
			Vec3f ambientComponent = computeAmbient(scene);
			r = ambientComponent.x;
			g = ambientComponent.y;
			b = ambientComponent.z;

			// Compute diffuse and specular components per light source
			for (int light = 0; light < numberOfLights; light++)
			{
				bool shadow = false;
				PointLight currentLight = scene.point_lights[light];

				float lightToCam = findDistance(currentLight.position, currentCamera.position);

				// Calculate light direction vector
				Vec3f wi = subtract(currentLight.position, this->intersectionPoint);

				wi = normalizeVector(wi);

				// Calculate wi_epsion
				Vec3f wi_epsilon;
				wi_epsilon.x = wi.x * scene.shadow_ray_epsilon;
				wi_epsilon.y = wi.y * scene.shadow_ray_epsilon;
				wi_epsilon.z = wi.z * scene.shadow_ray_epsilon;

				// Create shadowRay
				Vec3f shadow_e = add(this->intersectionPoint, wi_epsilon);
				Vec3f shadow_d = wi;
				Ray shadowRay(shadow_e, shadow_d);

				// float tLight = findLength(wi);
				float tLight = subtract(currentLight.position, shadowRay.e).x / shadowRay.d.x;

				Vec3f halfVector = subtract(wi, this->d);
				halfVector = normalizeVector(halfVector);

				for (int sphereNumber = 0; sphereNumber < numberOfSpheres; sphereNumber++)
				{
					Sphere currentSphere = scene.spheres[sphereNumber];
					Vec3f center = scene.vertex_data[currentSphere.center_vertex_id - 1];
					float radius = currentSphere.radius;

					shadowRay.sphereIntersection(currentSphere, sphereNumber, center);
				}
				for (int triangleNumber = 0; triangleNumber < numberOfTriangles; triangleNumber++)
				{
					Triangle currentTriangle = scene.triangles[triangleNumber];
					Vec3f v0 = scene.vertex_data[currentTriangle.indices.v0_id - 1];
					Vec3f v1 = scene.vertex_data[currentTriangle.indices.v1_id - 1];
					Vec3f v2 = scene.vertex_data[currentTriangle.indices.v2_id - 1];

					shadowRay.triangleIntersection(currentTriangle, v0, v1, v2, triangleNumber);
				}

				if (shadowRay.objID != -1 && tLight > shadowRay.minT)
				{
					shadow = true;
				}
				if (!shadow)
				{
					for (int meshNumber = 0; meshNumber < numberOfMeshes; meshNumber++)
					{
						Mesh currentMesh = scene.meshes[meshNumber];
						shadowRay.meshIntersection(currentMesh, scene, meshNumber);
					}
					if (shadowRay.objID != -1 && tLight > shadowRay.minT)
					{
						shadow = true;
					}
				}
				if (!shadow || (shadow && lightToCam == 0))
				{
					// no shadow, there will be contribution from the light source
					Vec3f irradiance = findIrradiance(currentLight);

					Vec3f diffuse = findDiffuse(currentLight, scene, irradiance);

					Vec3f specular = findSpecular(currentLight, scene, irradiance, halfVector);

					r += diffuse.x + specular.x;
					g += diffuse.y + specular.y;
					b += diffuse.z + specular.z;
				}
			}

			/************MIRRORNESS**********/

			bool mirrorness = isMirror(scene);

			Vec3f reflection;
			reflection.x = 0;
			reflection.y = 0;
			reflection.z = 0;

			if (maxDepth > 0 && mirrorness)
			{
				float wi = -2 * dot(this->d, this->surfaceNormal);
				Vec3f normal_wi;
				normal_wi.x = surfaceNormal.x * wi + d.x;
				normal_wi.y = surfaceNormal.y * wi + d.y;
				normal_wi.z = surfaceNormal.z * wi + d.z;

				normal_wi = normalizeVector(normal_wi);

				Vec3f wi_epsilon;
				wi_epsilon.x = normal_wi.x * scene.shadow_ray_epsilon;
				wi_epsilon.y = normal_wi.y * scene.shadow_ray_epsilon;
				wi_epsilon.z = normal_wi.z * scene.shadow_ray_epsilon;

				Vec3f e = add(this->intersectionPoint, wi_epsilon);
				Vec3f d = normal_wi;
				Ray reflectionRay(e, d);

				int numberOfSpheres = scene.spheres.size();
				int numberOfTriangles = scene.triangles.size();
				int numberOfMeshes = scene.meshes.size();

				/***********INTERSECT WITH SPHERES***********/
				for (int sphereNumber = 0; sphereNumber < numberOfSpheres; sphereNumber++)
				{
					Sphere currentSphere = scene.spheres[sphereNumber];
					Vec3f center = scene.vertex_data[currentSphere.center_vertex_id - 1];
					float radius = currentSphere.radius;
					reflectionRay.sphereIntersection(currentSphere, sphereNumber, center);
				}

				/***********INTERSECT WITH TRIANGLES***********/
				for (int triangleNumber = 0; triangleNumber < numberOfTriangles; triangleNumber++)
				{
					Triangle currentTriangle = scene.triangles[triangleNumber];
					Vec3f v0 = scene.vertex_data[currentTriangle.indices.v0_id - 1];
					Vec3f v1 = scene.vertex_data[currentTriangle.indices.v1_id - 1];
					Vec3f v2 = scene.vertex_data[currentTriangle.indices.v2_id - 1];
					reflectionRay.triangleIntersection(currentTriangle, v0, v1, v2, triangleNumber);
				}

				/**********INTERSECT WITH MESHES***********/
				for (int meshNumber = 0; meshNumber < numberOfMeshes; meshNumber++)
				{
					Mesh currentMesh = scene.meshes[meshNumber];
					reflectionRay.meshIntersection(currentMesh, scene, meshNumber);
				}

				if (!(reflectionRay.sphere && this->sphere && reflectionRay.objID == this->objID))
				{
					reflection = reflectionRay.computeColor(scene, currentCamera, (maxDepth - 1));
				}
				else if (!(reflectionRay.triangle && this->sphere && reflectionRay.objID == this->objID))
				{
					reflection = reflectionRay.computeColor(scene, currentCamera, (maxDepth - 1));
				}
				else if (!(reflectionRay.mesh && this->mesh && reflectionRay.objID == this->objID))
				{
					reflection = reflectionRay.computeColor(scene, currentCamera, (maxDepth - 1));
				}

				r += reflection.x * scene.materials[material_id - 1].mirror.x;
				g += reflection.y * scene.materials[material_id - 1].mirror.y;
				b += reflection.z * scene.materials[material_id - 1].mirror.z;
			}
		}
		else
		{
			if (maxDepth == scene.max_recursion_depth)
			{
				r = scene.background_color.x;
				g = scene.background_color.y;
				b = scene.background_color.z;
			}
			else
			{
				r = 0;
				g = 0;
				b = 0;
			}
		}

		if (r > 255)
		{
			r = 255;
		}
		else
		{
			r = round(r);
		}
		if (g > 255)
		{
			g = 255;
		}
		else
		{
			g = round(g);
		}
		if (b > 255)
		{
			b = 255;
		}
		else
		{
			b = round(b);
		}
		color.x = r;
		color.y = g;
		color.z = b;

		return color;
	}

	Vec3f findIrradiance(const PointLight &currentLight)
	{
		Vec3f irradiance;
		Vec3f d = subtract(currentLight.position, intersectionPoint);
		float d_square = dot(d, d);

		if (d_square != 0.0)
		{
			irradiance.x = currentLight.intensity.x / d_square;
			irradiance.y = currentLight.intensity.y / d_square;
			irradiance.z = currentLight.intensity.z / d_square;
		}
		return irradiance;
	}
	const Vec3f findDiffuse(const PointLight &currentLight, const Scene &scene, Vec3f irradiance)
	{
		Vec3f diffuse;

		Vec3f l = subtract(currentLight.position, intersectionPoint);
		l = normalizeVector(l);

		float dotPro = dot(l, this->surfaceNormal);
		if (dotPro < 0)
		{
			dotPro = 0;
		}

		diffuse.x = scene.materials[material_id - 1].diffuse.x * dotPro * irradiance.x;
		diffuse.y = scene.materials[material_id - 1].diffuse.y * dotPro * irradiance.y;
		diffuse.z = scene.materials[material_id - 1].diffuse.z * dotPro * irradiance.z;

		return diffuse;
	}

	Vec3f findSpecular(const PointLight &currentLight, const Scene &scene, Vec3f irradiance, Vec3f halfVector)
	{
		Vec3f specular;
		Material material = scene.materials[this->material_id - 1];

		float dotPro = dot(this->surfaceNormal, halfVector);
		if (dotPro < 0)
		{
			dotPro = 0;
		}

		specular.x = material.specular.x * pow(dotPro, material.phong_exponent) * irradiance.x;
		specular.y = material.specular.y * pow(dotPro, material.phong_exponent) * irradiance.y;
		specular.z = material.specular.z * pow(dotPro, material.phong_exponent) * irradiance.z;

		return specular;
	}
	bool isMirror(const Scene &scene)
	{
		if (scene.materials[this->material_id - 1].mirror.x > 0 || scene.materials[this->material_id - 1].mirror.y > 0 || scene.materials[this->material_id - 1].mirror.z > 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

int main(int argc, char *argv[])
{
	// Sample usage for reading an XML scene file
	Scene scene;

	scene.loadFromXml(argv[1]);

	int numberOfCameras = scene.cameras.size();

	/**********FOR EACH PIXEL IN EACH CAMERA*************/
	for (int cameraNumber = 0; cameraNumber < numberOfCameras; cameraNumber++)
	{
		Camera currentCamera = scene.cameras[cameraNumber];
		int width = currentCamera.image_width;
		int height = currentCamera.image_height;

		unsigned char *image = new unsigned char[width * height * 3];
		int pixelNumber = 0;

		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				/*********GENERATE RAYS*********/
				Ray ray(currentCamera, i, j);

				int numberOfSpheres = scene.spheres.size();
				int numberOfTriangles = scene.triangles.size();
				int numberOfMeshes = scene.meshes.size();

				/***********INTERSECT WITH SPHERES***********/
				for (int sphereNumber = 0; sphereNumber < numberOfSpheres; sphereNumber++)
				{
					Sphere currentSphere = scene.spheres[sphereNumber];
					Vec3f center = scene.vertex_data[currentSphere.center_vertex_id - 1];

					ray.sphereIntersection(currentSphere, sphereNumber, center);
				}

				/***********INTERSECT WITH TRIANGLES***********/
				for (int triangleNumber = 0; triangleNumber < numberOfTriangles; triangleNumber++)
				{
					Triangle currentTriangle = scene.triangles[triangleNumber];
					Vec3f v0 = scene.vertex_data[currentTriangle.indices.v0_id - 1];
					Vec3f v1 = scene.vertex_data[currentTriangle.indices.v1_id - 1];
					Vec3f v2 = scene.vertex_data[currentTriangle.indices.v2_id - 1];

					ray.triangleIntersection(currentTriangle, v0, v1, v2, triangleNumber);
				}

				/**********INTERSECT WITH MESHES***********/
				for (int meshNumber = 0; meshNumber < numberOfMeshes; meshNumber++)
				{
					Mesh currentMesh = scene.meshes[meshNumber];
					ray.meshIntersection(currentMesh, scene, meshNumber);
				}

				Vec3f pixelColor = ray.computeColor(scene, currentCamera, (scene.max_recursion_depth));

				image[pixelNumber] = pixelColor.x;
				image[pixelNumber + 1] = pixelColor.y;
				image[pixelNumber + 2] = pixelColor.z;
				pixelNumber += 3;
			}
		}

		write_ppm(currentCamera.image_name.c_str(), image, width, height);
	}
}