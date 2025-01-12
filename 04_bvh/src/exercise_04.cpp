#include <cglib/rt/bvh.h>

#include <cglib/rt/triangle_soup.h>

#include <cglib/core/image.h>
#include <complex>
#include <map>

/*
 * Create a 1 dimensional normalized gauss kernel
 *
 * Parameters:
 *  - sigma:       the parameter of the gaussian
 *  - kernel_size: the size of the kernel (has to be odd)
 *  - kernel:      an array with size kernel_size elements that
 *                 has to be filled with the kernel values
 *
 */
void Image::create_gaussian_kernel_1d(
    float sigma,
    int kernel_size,
    float* kernel) {
    cg_assert(kernel_size%2==1);

    // TODO: calculate filter values as described on the exercise sheet.
    // Make sure your kernel is normalized
    float sum = 0;
    float normalization = 1 / (sigma * glm::sqrt(2 * std::numbers::pi));
    int mid_point = int(kernel_size / 2);

    for (int i = 0; i < kernel_size; ++i) {
        kernel[i] = normalization * glm::exp(-glm::pow(i - mid_point, 2) / (2 * sigma * sigma));
        sum = sum + kernel[i];
    }
    // normalizing
    for (int i = 0; i < kernel_size; ++i) {
        kernel[i] = kernel[i] / sum;
    }
}

/*
 * Create a 2 dimensional quadratic and normalized gauss kernel
 *
 * Parameters:
 *  - sigma:       the parameter of the gaussian
 *  - kernel_size: the size of the kernel (has to be odd)
 *  - kernel:      an array with kernel_size*kernel_size elements that
 *                 has to be filled with the kernel values
 */
void Image::create_gaussian_kernel_2d(
    float sigma,
    int kernel_size,
    float* kernel) {
    cg_assert(kernel_size%2==1);

    // TODO: calculate filter values as described on the exercise sheet.
    // Make sure your kernel is normalized
    float normalization = 1 / (sigma * sigma * 2 * std::numbers::pi);
    int mid_point = int(kernel_size / 2);
    float sum = 0.f;
    for (int j = 0; j < kernel_size; ++j) {
        for (int i = 0; i < kernel_size; ++i) {
            kernel[i + j * kernel_size] = normalization * glm::exp(
                -(glm::pow(i - mid_point, 2) + glm::pow(j - mid_point, 2)) / (2 * sigma * sigma));
            sum = sum + kernel[i + j * kernel_size];
        }
    }
    // normalization
    for (int j = 0; j < kernel_size; ++j) {
        for (int i = 0; i < kernel_size; ++i) {
            kernel[i + j * kernel_size] = kernel[i + j * kernel_size] / sum;
        }
    }
}

/*
 * Convolve an image with a 2d filter kernel
 *
 * Parameters:
 *  - kernel_size: the size of the 2d-kernel
 *  - kernel:      the 2d-kernel with kernel_size*kernel_size elements
 *  - wrap_mode:   needs to be known to handle repeating
 *                 textures correctly
 */
void Image::filter(Image* target, int kernel_size, float* kernel, WrapMode wrap_mode) const {
    cg_assert(kernel_size%2==1 && "kernel size should be odd.");
    cg_assert(kernel_size > 0 && "kernel size should be greater than 0.");
    cg_assert(target);
    cg_assert(target->getWidth() == m_width && target->getHeight() == m_height);

    int midpoint = int(kernel_size / 2);

    for (int i = 0; i < target->getWidth(); ++i) {
        for (int j = 0; j < target->getHeight(); ++j) {
            glm::vec4 pixel = glm::vec4(0.f);
            for (int k = 0; k < kernel_size; ++k) {
                for (int l = 0; l < kernel_size; ++l) {
                    pixel += this->getPixel(i - k + midpoint, j - l + midpoint, wrap_mode) * kernel[l + k *
                        kernel_size];
                }
            }
            target->setPixel(i, j, pixel);
        }
    }
}

/*
 * Convolve an image with a separable 1d filter kernel
 *
 * Parameters:
 *  - kernel_size: the size of the 1d kernel
 *  - kernel:      the 1d-kernel with kernel_size elements
 *  - wrap_mode:   needs to be known to handle repeating
 *                 textures correctly
 */
void Image::filter_separable(Image* target, int kernel_size, float* kernel, WrapMode wrap_mode) const {
    cg_assert(kernel_size%2==1 && "kernel size should be odd.");
    cg_assert(kernel_size > 0 && "kernel size should be greater than 0.");
    cg_assert(target);
    cg_assert(target->getWidth() == m_width && target->getHeight() == m_height);

    // TODO: realize the 2d convolution with two
    // convolutions of the image with a 1d-kernel.
    // convolve the image horizontally and then convolve
    // the result vertically (or vise-versa).
    //
    // use the methods getPixel(x, y, wrap_mode) and
    // setPixel(x, y, value) to get and set pixels of an image

    int midpoint = int(kernel_size / 2);
    for (int i = 0; i < target->getWidth(); ++i) {
        for (int j = 0; j < target->getHeight(); ++j) {
            glm::vec4 pixel = glm::vec4(0.f);
            for (int k = 0; k < kernel_size; ++k) {
                pixel += this->getPixel(i - k + midpoint, j, wrap_mode) * kernel[k];
            }
            target->setPixel(i, j, pixel);
        }
    }
    for (int i = 0; i < target->getWidth(); ++i) {
        for (int j = 0; j < target->getHeight(); ++j) {
            glm::vec4 pixel = glm::vec4(0.f);
            for (int k = 0; k < kernel_size; ++k) {
                pixel += target->getPixel(i, j - k + midpoint, wrap_mode) * kernel[k];
            }
            target->setPixel(i, j, pixel);
        }
    }
}

/**
 * Reorder triangle indices in the vector triangle_indices
 * in the range [first_triangle_idx, first_triangle_idx+num_triangles-1]
 * so that the range is split in two sets where all triangles in the first set
 * are "less than equal" than the median, and all triangles in the second set
 * are "greater than equal" the median.
 *
 * Ordering ("less than") is defined by the ordering of triangle
 * bounding box centers along the given axis.
 *
 * Triangle indices within a set need not be sorted.
 *
 * The resulting sets must have an equal number of elements if num_triangles
 * is even. Otherwise, one of the sets must have one more element.
 *
 * For example, 8 triangles must be split 4-4. 7 Triangles must be split
 * 4-3 or 3-4.
 *
 * Parameters:
 *  - first_triangle_idx: The index of the first triangle in the given range.
 *  - num_triangles:      The number of triangles in the range.
 *  - axis:               The sort axis. 0 is X, 1 is Y, 2 is Z.
 *
 * Return value:
 *  - The number of triangles in the first set.
 */
int BVH::reorder_triangles_median(
    int first_triangle_idx,
    int num_triangles,
    int axis) {
    cg_assert(first_triangle_idx < static_cast<int>(triangle_indices.size()));
    cg_assert(first_triangle_idx >= 0);
    cg_assert(num_triangles <= static_cast<int>(triangle_indices.size() - first_triangle_idx));
    cg_assert(num_triangles > 1);
    cg_assert(axis >= 0);
    cg_assert(axis < 3);

    std::vector<std::pair<int, float>> midPoints_along_axis;
    // as we need the middle of each aabb we need to first creat those aabb containers
    for (int i = first_triangle_idx; i < num_triangles + first_triangle_idx; ++i) {
        AABB aabb;
        aabb.extend(this->triangle_soup.vertices[this->triangle_indices[i] * 3]);
        aabb.extend(this->triangle_soup.vertices[this->triangle_indices[i] * 3 + 1]);
        aabb.extend(this->triangle_soup.vertices[this->triangle_indices[i] * 3 + 2]);
        glm::vec3 diff = aabb.max - aabb.min;
        diff /= 2;
        glm::vec3 end = aabb.min + diff;
        midPoints_along_axis.emplace_back(this->triangle_indices[i], end[axis]);
    }
    //sort
    std::sort(midPoints_along_axis.begin(), midPoints_along_axis.end(),
              [](const std::pair<int, float>& a, const std::pair<int, float>& b) -> bool {
                  return a.second > b.second;
              });

    for (long unsigned int i = 0; i < midPoints_along_axis.size(); ++i) {
        this->triangle_indices[i + first_triangle_idx] = midPoints_along_axis[i].first;
    }
    return static_cast<int>(midPoints_along_axis.size() / 2);
}

/*
 * Build a BVH recursively using the object median split heuristic.
 *
 * This method must first fully initialize the current node, and then
 * potentially split it.
 *
 * A node must not be split if it contains MAX_TRIANGLES_IN_LEAF triangles or
 * less. No leaf node may be empty. All nodes must have either two or no
 * children.
 *
 * Use reorder_triangles_median to perform the split in triangle_indices.
 * Split along X for depth 0. Then, proceed in the order Y, Z, X, Y, Z, X, Y, ..
 *
 * Parameters:
 *  - node_idx:           The index of the node to be split.
 *  - first_triangle_idx: An index into the array triangle_indices. It points
 *                        to the first triangle contained in the current node.
 *  - num_triangles:      The number of triangles contained in the current node.
 */
void BVH::
build_bvh(int node_idx, int first_triangle_idx, int num_triangles, int depth) {
    cg_assert(num_triangles > 0);
    cg_assert(node_idx >= 0);
    cg_assert(node_idx < static_cast<int>(nodes.size()));
    cg_assert(depth >= 0);

    nodes[node_idx].triangle_idx = first_triangle_idx;
    nodes[node_idx].num_triangles = num_triangles;
    nodes[node_idx].aabb.min = glm::vec3(FLT_MAX);
    nodes[node_idx].aabb.max = glm::vec3(-FLT_MAX);
    nodes[node_idx].left = -1;
    nodes[node_idx].right = -1;
    cg_assert(num_triangles>0);

    if (num_triangles <= MAX_TRIANGLES_IN_LEAF) {
        for (int i = first_triangle_idx; i < num_triangles + first_triangle_idx; ++i) {
            nodes[node_idx].aabb.extend(triangle_soup.vertices[triangle_indices[i] * 3]);
            nodes[node_idx].aabb.extend(triangle_soup.vertices[triangle_indices[i] * 3 + 1]);
            nodes[node_idx].aabb.extend(triangle_soup.vertices[triangle_indices[i] * 3 + 2]);
        }
        return;
    }

    nodes.push_back(Node{});
    nodes.push_back(Node{});
    nodes[node_idx].left = nodes.size() - 2;
    nodes[node_idx].right = nodes.size() - 1;
    int left_first_idx = this->reorder_triangles_median(first_triangle_idx, num_triangles, depth % 3);
    this->build_bvh(nodes[node_idx].left, first_triangle_idx, left_first_idx, depth + 1);
    this->build_bvh(nodes[node_idx].right, first_triangle_idx + left_first_idx, num_triangles - left_first_idx,
                    depth + 1);

    nodes[node_idx].aabb.extend(nodes[nodes[node_idx].right].aabb.min);
    nodes[node_idx].aabb.extend(nodes[nodes[node_idx].right].aabb.max);
    nodes[node_idx].aabb.extend(nodes[nodes[node_idx].left].aabb.min);
    nodes[node_idx].aabb.extend(nodes[nodes[node_idx].left].aabb.max);
}

/*
 * Intersect the BVH recursively, returning the nearest intersection if
 * there is one.
 *
 * Caution: BVH nodes may overlap.
 *
 * Parameters:
 *  - ray:                  The ray to intersect the BVH with.
 *  - idx:                  The node to be intersected.
 *  - nearest_intersection: The distance to the intersection point, if an
 *                          intersection was found. Must not be changed
 *                          otherwise.
 *  - isect:                The intersection, if one was found. Must not be
 *                          changed otherwise.
 *
 * Return value:
 *  true if an intersection was found, false otherwise.
 */
bool BVH::
intersect_recursive(const Ray& ray, int idx, float* nearest_intersection, Intersection* isect) const {
    cg_assert(nearest_intersection);
    cg_assert(isect);
    cg_assert(idx >= 0);
    cg_assert(idx < static_cast<int>(nodes.size()));

    const Node& n = nodes[idx];

    // This is a leaf node. Intersect all triangles.
    if (n.left < 0) {
        glm::vec3 bary(0.f);
        bool hit = false;
        for (int i = 0; i < n.num_triangles; i++) {
            int x = triangle_indices[n.triangle_idx + i];
            float dist;
            glm::vec3 b;
            if (intersect_triangle(ray.origin, ray.direction,
                                   triangle_soup.vertices[x * 3 + 0],
                                   triangle_soup.vertices[x * 3 + 1],
                                   triangle_soup.vertices[x * 3 + 2],
                                   b, dist)) {
                hit = true;
                if (dist <= *nearest_intersection) {
                    *nearest_intersection = dist;
                    bary = b;
                    cg_assert(x >= 0);
                    if (isect)
                        triangle_soup.fill_intersection(isect, x, *nearest_intersection, bary);
                }
            }
        }
        return hit;
    }

    // This is an inner node. Recurse into child nodes.
    else {
        auto tmin_l = -FLT_MAX;
        auto tmin_r = -FLT_MAX;
        auto tmax_l = FLT_MAX;
        auto tmax_r = FLT_MAX;

        bool isect_l = nodes[n.left].aabb.intersect(ray, tmin_l, tmax_l);
        bool isect_r = nodes[n.right].aabb.intersect(ray, tmin_r, tmax_r);
        if (isect_l || isect_r) {
            // if both have been hit
            if (isect_l && isect_r) {
                bool hit = false;
                if (tmin_l <= tmin_r) {
                    hit = intersect_recursive(ray, n.left, nearest_intersection, isect);
                }
                else {
                    hit = intersect_recursive(ray, n.right, nearest_intersection, isect);
                }
                // if the nearer hasnt got an result
                if (!hit) {
                    if (tmin_l < tmin_r) {
                        hit = intersect_recursive(ray, n.right, nearest_intersection, isect);
                    }
                    else {
                        hit = intersect_recursive(ray, n.left, nearest_intersection, isect);
                    }
                }
                return hit;
            }
            // if only one has been hit
            if (isect_l) {
                return intersect_recursive(ray, n.left, nearest_intersection, isect);
            }
            return intersect_recursive(ray, n.right, nearest_intersection, isect);
        }
    }

    return false;
}


// CG_REVISION 8c58412a25ac2367c053bfa840ee81b320bdd315
