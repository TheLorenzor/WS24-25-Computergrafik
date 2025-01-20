#include <iostream>
#include <cglib/core/image.h>

#include <cglib/gl/scene_graph.h>
#include <cglib/gl/prefilter_envmap.h>

#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtx/io.hpp>

/*
 * Create a sphere flake node and its 5 children. 
 * Recurse if number_of_remaining_recursions > 0.
 */
std::shared_ptr<SceneGraphNode> buildSphereFlakeSceneGraph(
    const std::shared_ptr<GLObjModel>& model, // the geometric model associated with this node (a sphere)
    float sphere_model_radius,
    // the radius of the sphere in local coordinates (hint: should remain constant throughout the recursion)
    float size_factor, // growth/shrink factor of the spheres for each recursion level
    int number_of_remaining_recursions) // the number of remaining recursions
{
    // TODO: create root SceneGraphNode and set model
    auto root = std::make_shared<SceneGraphNode>();
    root->model = model;

    // TODO: only create children if there are remaining recursions
    if (number_of_remaining_recursions <= 0) {
        return root;
    }
    std::vector<glm::vec3> directions = {
        {1, 0, 0}, // +X
        {-1, 0, 0}, // -X
        {0, 0, 1}, // +z
        {0, 0, -1}, // -z
        {0, 1, 0}, // +y
    };
    // TODO: create 5 child nodes
    number_of_remaining_recursions--;
    for (int i = 0; i < 5; ++i) {
        // TODO: create sphere flake sub graph for each child
        auto subgraph = buildSphereFlakeSceneGraph(model,
                                                   sphere_model_radius, size_factor,
                                                   number_of_remaining_recursions);

        // TODO: compute transformation matrix from the child node to the parent node (this function call's root node)
        glm::mat4 translation = glm::translate(directions[i]*(sphere_model_radius+size_factor));
        glm::mat4 scale = glm::scale(glm::vec3(size_factor));
        glm::vec3 rot_direction = glm::vec3(directions[i].z,directions[i].y,directions[i].x);
        glm::mat4 rotate = glm::rotate(glm::half_pi<float>(),rot_direction);
        subgraph->node_to_parent =translation*rotate*scale;
        subgraph->parent_to_node = glm::inverse(subgraph->node_to_parent);

        // TODO: add each subgraph to the children of this function call's root node
        root->children.push_back(subgraph);
    }

    return root;
}


/*
 * Traverse the scene graph and collect all models.
 * Recursivly computes the world transformation for each model.
 */
void SceneGraphNode::collectTransformedModels(
    std::vector<TransformedModel>& transformed_models, // list that contains all models when traversal is complete
    const glm::mat4& parent_to_world, // the world transformation of the parent of this node
    const glm::mat4& world_to_parent) // the inverse world transformation of this node's parent
const {
    // TODO: compute node_to_world and world_to_node transformation matrices for this node
    glm::mat4 node_to_world =  this->node_to_parent* parent_to_world;
    glm::mat4 world_to_node = world_to_parent * this->parent_to_node;
    //std::cout <<this->node_to_parent << std::endl;
    //std::cout << "---------------------------------------" << std::endl;

    // TODO: add this node's model to the list of transformed models
    transformed_models.push_back(    TransformedModel(node_to_world, world_to_node,this->model.get()));
    // TODO: recursively transform and add the models of all children (subgraphs)
    for (auto& child : this->children) {
        child->collectTransformedModels(transformed_models, world_to_node, world_to_node);
    }
}


/*
 * Perform animation by traversing the scene graph recursively
 * and rotating each scene graph node around the y-axis (in object space).
 */
void animateSphereFlake(
    SceneGraphNode& node, // the current node of the traversal
    float angle_increment) // the incremental rotation angle
{
    // TODO: compute incremental rotation matrix for the given node node
    glm::mat4 rotation = glm::rotate(angle_increment, glm::vec3(0,1,0));

    // TODO: compute the parent-relative transformation matrices for this node
    node.node_to_parent = rotation*node.node_to_parent;
    node.parent_to_node = glm::transpose(node.node_to_parent);

    // TODO: recursively animate child subgraphs
    for (auto& child : node.children) {
        animateSphereFlake(*child,angle_increment);
    }
}


/*
 * compute a prefiltered environment map for diffuse
 * illumination for the given environment map
 */
std::shared_ptr<Image>
prefilter_environment_diffuse(Image const& img) {
    const int width = img.getWidth();
    const int height = img.getHeight();
    auto filtered = std::make_shared<Image>(width, height);

    // For all texels in the envmap...
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // TODO: compute normal direction

            // ... integrate over all incident directions.
            for (int dy = 0; dy < height; ++dy) {
                for (int dx = 0; dx < width; ++dx) {
                    // TODO: compute incident direction
                    // TODO: accumulate samples
                }
            }

            // TODO: write filtered value
            //filtered->setPixel(...);
        }
    }

    return filtered;
}

/*
 * compute a prefiltered environment map for specular
 * illumination for the given environment map and
 * phong exponent
 */
std::shared_ptr<Image>
prefilter_environment_specular(Image const& img, float n) {
    const int width = img.getWidth();
    const int height = img.getHeight();

    auto filtered = std::make_shared<Image>(width, height);

    // For all texels in the envmap...
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // TODO: compute reflection direction

            // ... integrate over all incident directions.
            for (int dy = 0; dy < height; ++dy) {
                for (int dx = 0; dx < width; ++dx) {
                    // TODO: compute incident direction

                    // TODO: accumulate samples
                }
            }

            // TODO: write filtered value
            //filtered->setPixel(...);
        }
    }

    return filtered;
}

// CG_REVISION d4ab32bd208749f2d2b1439e25d16e642b039298
