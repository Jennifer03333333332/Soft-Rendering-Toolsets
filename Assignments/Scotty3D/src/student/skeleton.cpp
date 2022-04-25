
#include "../scene/skeleton.h"
//bind position:joint space -> skeleton space
Mat4 Joint::joint_to_bind() const {

    // TODO(Animation): Task 2

    // Return a matrix transforming points in the joint space of this joint 
    // to points in skeleton space in bind position.

    // Bind position implies that all joints have pose = Vec3{0.0f}

    Mat4 iter = Mat4::I;
    //? should I consider the current extent?

    Joint* cur = parent;
    while(cur){// !cur->is_root()
        //translate origin to the extent of the parents
        iter = Mat4::translate(cur->extent)*iter;
        cur = cur->parent;
    }
    return iter;
}

//pose: joint space -> skeleton space
Mat4 Joint::joint_to_posed() const {

    // TODO(Animation): Task 2

    // Return a matrix transforming points in the joint space of this joint
    // to points in skeleton space in posed position.

    // Posed position implies that you should take into account the joint 
    // poses along with the extents.

    //why we don't need translate for the current joint?hmmm...
    Mat4 iter = Mat4::euler(pose);

    //error 1 Mat4 iter = Mat4::euler(pose)*Mat4::translate(extent);


    Joint* cur = parent;
    while(cur){
        //question: rotate first or translate first? because rotate must around the origin, so translate first?
        iter = Mat4::translate(cur->extent)*iter;
        iter = Mat4::euler(cur->pose)*iter;//R:=Rz⋅Ry⋅Rx
        
        cur = cur->parent;
    }
    return iter;
}
//Vec3,endpoint's bind: joint -> world space
Vec3 Skeleton::end_of(Joint* j) {

    // TODO(Animation): Task 2

    // Return the bind position of the endpoint of joint j in world space.
    // This should take into account Skeleton::base_pos.
    Vec3 world_bindpos = base_pos;
    //what's the endpoint of joint?
    Vec3 endpoint = j->extent;
    Vec3 skeleton_bindpos = j->joint_to_bind()*endpoint;

    //translate equals to add though
    //Mat4 skeleton2world = Mat4::translate(base_pos);
    //world_bindpos = skeleton2world*skeleton_bindpos;

    world_bindpos += skeleton_bindpos;
    return world_bindpos;
}

//Vec3, pose: joint -> world space
Vec3 Skeleton::posed_end_of(Joint* j) {

    // TODO(Animation): Task 2

    // Return the posed position of the endpoint of joint j in world space.
    // This should take into account Skeleton::base_pos.
    Vec3 world_posepos = base_pos;

    Vec3 endpoint = j->extent;
    Vec3 skeleton_posepos = j->joint_to_posed()*endpoint;

    //translate equals to add though
    // Mat4 skeleton2world = Mat4::translate(base_pos);
    // world_posepos = skeleton2world*skeleton_posepos;
    world_posepos += skeleton_posepos;
    return world_posepos;
}

//Matrix: bindpos from joint space ->world space
Mat4 Skeleton::joint_to_bind(const Joint* j) const {

    // TODO(Animation): Task 2

    // Return a matrix transforming points in joint j's space to world space in
    // bind position. This should take into account Skeleton::base_pos.
    
    // Calculate skeleton2world matrix
    Mat4 skeleton2world = Mat4::I;
    skeleton2world = skeleton2world.translate(base_pos);

    return skeleton2world*j->joint_to_bind();
}
//Matrix: pose from joint space ->world space
Mat4 Skeleton::joint_to_posed(const Joint* j) const {

    // TODO(Animation): Task 2

    // Return a matrix transforming points in joint j's space to world space with
    // poses. This should take into account Skeleton::base_pos.
    Mat4 skeleton2world = Mat4::I;
    skeleton2world = skeleton2world.translate(base_pos);
    return skeleton2world*j->joint_to_posed();
}

void Joint::compute_gradient(Vec3 target, Vec3 current) {

    // TODO(Animation): Task 2

    // Computes the gradient of IK energy for this joint and, should be called
    // recursively upward in the heirarchy. Each call should storing the result
    // in the angle_gradient for this joint.

    // 1 jacobian of θ = r x p
    // r is the axis of rotation in the current joint space.
    Vec3 x_axis = joint_to_posed()*Vec3(1,0,0);//joint_to_posed()*??
    Vec3 y_axis = joint_to_posed()*Vec3(0,1,0);
    Vec3 z_axis = joint_to_posed()*Vec3(0,0,1);

    // p is the vector from the base of joint i to the end point of the target joint.
    // joint_to_posed() * Vec3(0,0,0) is the base of joint i
    Vec3 p = current - joint_to_posed() * Vec3(0,0,0);//in skeleton space //error ! current - joint_to_posed() * Vec3(0,0,0);
    Vec3 Jacobian_x = cross(x_axis,p);//for this joint, this rotation axis
    Vec3 Jacobian_y = cross(y_axis,p);
    Vec3 Jacobian_z = cross(z_axis,p);


    angle_gradient.x += dot(Jacobian_x,current - target);//delta_f
    angle_gradient.y += dot(Jacobian_y,current - target);
    angle_gradient.z += dot(Jacobian_z,current - target);

    //how to find all my parents?
    Joint* cur = parent;
    while(cur){
        x_axis = cur->joint_to_posed()*Vec3(1,0,0);//joint_to_posed()*??
        y_axis = cur->joint_to_posed()*Vec3(0,1,0);
        z_axis = cur->joint_to_posed()*Vec3(0,0,1);
        // p is the vector from the base of joint i to the end point of the target joint.
        //? target or current
        p = current - cur->joint_to_posed() * Vec3(0,0,0);//in skeleton space 
        Jacobian_x = cross(x_axis,p);//for this joint, this rotation axis
        Jacobian_y = cross(y_axis,p);
        Jacobian_z = cross(z_axis,p);

        cur->angle_gradient.x += dot(Jacobian_x,current - target);//delta_f
        cur->angle_gradient.y += dot(Jacobian_y,current - target);
        cur->angle_gradient.z += dot(Jacobian_z,current - target);
        cur = cur->parent;
    }
    // q: Target is the target position of the IK handle in skeleton space.
    // p(theta t): Current is the end position of the IK'd joint in skeleton space.
}

void Skeleton::step_ik(std::vector<IK_Handle*> active_handles) {

    // TODO(Animation): Task 2
    
    //float cost_func = 0;

    //stop after 
    float frames = 50.0f;
    float tau = 1 / frames;//timestep
    //where is the alpha?
    while(frames--){
        //one joint could be used in many IK_Handle
        for(size_t i = 0; i<active_handles.size();i++){
            Joint* cur_joint = active_handles[i]->joint;
            Vec3 end_pos_thisjoint = cur_joint->joint_to_posed()*cur_joint->extent;//should in skeleton space
            cur_joint->compute_gradient(active_handles[i]->target, end_pos_thisjoint);
        }
        //for_joints, add a little to pose
        for_joints([&](Joint* cur_joint) {
            //get the angle_gradient
            cur_joint->pose = cur_joint->pose - tau*cur_joint->angle_gradient;//Mat4::euler(cur_joint->angle_gradient)*cur_joint->pose
            //finally clean the angle_gradient
            cur_joint->angle_gradient = Vec3();
        });
    }
    // Do several iterations of Jacobian Transpose gradient descent for IK
}


//All in joint space? Depends
Vec3 closest_on_line_segment(Vec3 start, Vec3 end, Vec3 point) {

    // TODO(Animation): Task 3

    // Return the closest point to 'point' on the line segment from start to end

    // project to this line, see where's the point
    // if out, return the start/end

    Vec3 start_p = point - start;
    Vec3 start_end = end - start;
    //Situation 1: < start
    if(dot(start_p,start_end) <= 0)return start;

    //proj = Length * direction
    Vec3 proj = dot(start_p,start_end) / start_end.norm() * start_end / start_end.norm();//From start -> project point

    //Situation 2: > end
    if(proj.norm_squared() > start_end.norm_squared())return end;


    return proj;
}

void Skeleton::find_joints(const GL::Mesh& mesh, std::vector<std::vector<Joint*>>& map) {

    // TODO(Animation): Task 3

    // Construct a mapping: vertex index -> list of joints that should effect the vertex.
    // A joint should effect a vertex if it is within Joint::radius distance of the bone's line segment in bind position.
    // ? the verts is in skeleton space or world space?try world space
    const std::vector<GL::Mesh::Vert>& verts = mesh.verts();
    map.resize(verts.size());

    // For each i in [0, verts.size()), map[i] should contain the list of joints that
    // effect vertex i. Note that i is NOT Vert::id! i is the index in verts.

    for_joints([&](Joint* j) {
        // Judge the vertex's position - projection point's position
        for(size_t i = 0; i< verts.size(); i++){
            //How to find j's skeleton: nvm
            //bind world space -> bind joint space
            Vec3 vpos_in_jointspace = Mat4::inverse(joint_to_bind(j))  * (verts[i].pos);
            Vec3 proj = closest_on_line_segment(Vec3{0}, Vec3{0} + j->extent, vpos_in_jointspace);
            //joint j should effect this v:
            if((verts[i].pos - proj).norm() <= j->radius){
                map[i].push_back(j);
            }
        }
        

    });
}

void Skeleton::skin(const GL::Mesh& input, GL::Mesh& output,
                    const std::vector<std::vector<Joint*>>& map) {

    // TODO(Animation): Task 3

    // Apply bone poses & weights to the vertices of the input (bind position) mesh
    // and store the result in the output mesh. 
    // map was computed by find_joints, hence gives a mapping from vertex index to
    // the list of bones the vertex should be effected by.

    // Currently, this just copies the input to the output without modification.

    std::vector<GL::Mesh::Vert> verts = input.verts();

    for(size_t i = 0; i < verts.size(); i++) {

        // Skin vertex i. Note that its position is given in object bind space.
        float sum_of_inv_dis = 0;
        std::vector<float> inv_dist;
        //Try to save the performance
        //inv_dist.resize(map[i].size());//joint's quantity
        //for(auto j:map[i]){
        for(size_t j = 0; j < map[i].size(); j++){
            
            //Try All in world space
            //float distance_ij = (closest_on_line_segment(joint_to_bind(map[i][j])*Vec3{0}, end_of(map[i][j]), verts[i].pos) - verts[i].pos).norm();
            //Try All in joint space
           
            Mat4 world2jointspace_bind = Mat4::inverse(joint_to_bind(map[i][j]));
            Vec3 vpos_inJoint = world2jointspace_bind*verts[i].pos;
            float distance_ij = (closest_on_line_segment(Vec3{0}, Vec3{0} + map[i][j]->extent, vpos_inJoint) - vpos_inJoint).norm();
            //what if distance_ij == 0? the vertex's is exactly on this joint
            inv_dist.emplace_back(1.0f/distance_ij);//inv_dist[j] = 1.0f/distance_ij;
            sum_of_inv_dis += inv_dist[j];
        }
        if(inv_dist.empty()) continue;
        Vec3 sum_of_pos{0};
        for(size_t j = 0; j < map[i].size(); j++){
            Mat4 world2jointspace_bind = Mat4::inverse(joint_to_bind(map[i][j]));
            Vec3 v_ij = world2jointspace_bind*verts[i].pos;
            sum_of_pos += inv_dist[j]/sum_of_inv_dis * v_ij;// w_ij * v_ij

        }
        verts[i].pos = sum_of_pos;
    }

    std::vector<GL::Mesh::Index> idxs = input.indices();
    output.recreate(std::move(verts), std::move(idxs));
}
