#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {
//Rope has two members, vector<Mass *>masses and vector<Spring *>springs
    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.

        //Comment-in this part when you implement the constructor
        if(num_nodes==1)
            masses.push_back((new Mass(start,node_mass,false)));
        else
        {
            for (int i = 0; i < num_nodes; i++)
                masses.push_back(new Mass(start + i * (end - start) / (num_nodes - 1), node_mass, false));
        }
        for (auto &i : pinned_nodes)
            masses[i]->pinned = true;
        for(int i=0;i<num_nodes-1;i++)
            springs.push_back(new Spring(masses[i],masses[i+1],k));
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            const Vector2D diff=s->m2->position-s->m1->position;
            const Vector2D f=s->k*(diff.norm()-s->rest_length)*diff.unit();//f points to s->m2
            s->m1->forces+=f;
            s->m2->forces-=f;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces+=gravity*m->mass;
                m->forces-=0.01*m->velocity;
                const Vector2D a=m->forces/m->mass;
                m->velocity+=a*delta_t;
                m->position+=m->velocity*delta_t;
                // TODO (Part 2): Add global damping
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            const Vector2D diff=s->m2->position-s->m1->position;
            const Vector2D f=s->k*(diff.norm()-s->rest_length)*diff.unit();//f points to s->m2
            s->m1->forces+=f;
            s->m2->forces-=f;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass
                m->forces += gravity * m->mass;
                const Vector2D a = m->forces / m->mass;
                m->position += 0.99995 * (temp_position - m->last_position) + a * delta_t * delta_t;
                m->last_position = temp_position;
                // TODO (Part 4): Add global Verlet damping
            }
            m->forces=Vector2D(0.f,0.f);
        }
    }
}
