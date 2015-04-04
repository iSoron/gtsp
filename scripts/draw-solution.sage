#!/usr/bin/sage
from sage.plot.colors import red, white, blue
import sys

FIGURE_SIZE = 10
POINT_SIZE = FIGURE_SIZE

def turn_direction(a,b,c):
    x1,x2,x3=a[0],b[0],c[0]
    y1,y2,y3=a[1],b[1],c[1]
    return (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

def angle_distance(dx,dy):
    dx,dy = n(dx), n(dy)
    length = sqrt(dx*dx+dy*dy)
    return (dx / length, -length)

def chull_2d(points):
    """Computes the convex hull of a set of 2-dimensional points."""

    points = copy(points)

    # Find the point with lower y-coordinate
    pivot = min(points, key=(lambda x: (x[1],x[0])))
    points.remove(pivot)

    # Sorts all the points according to the angle
    points.sort(key=(lambda p : angle_distance(p[0]-pivot[0],p[1]-pivot[1])), reverse=True)

    # Add the two first points to the stack
    stack = [pivot,points[0]]

    k = 1
    while k < len(points):
        p0 = stack[len(stack)-2]
        p1 = stack[len(stack)-1]
        p2 = points[k]

        t = turn_direction(p0,p1,p2)
        
        # Left turn or straight
        if(t <= 0):
            stack.pop()

        # Right turn or straight
        if(t >= 0):
            stack.append(p2)
            k = k+1

    return stack

def read_problem_file(problem_file):
    (node_count, cluster_count, edges_count) = [ int(x) for x in problem_file.readline().split(' ') ]

    points = []
    all_points = []
    for i in range(cluster_count):
        points.append([])

    for i in range(node_count):
        (x,y,cluster) = [ int(float(x)) for x in problem_file.readline().split(' ') ]
        points[cluster].append((x,y))
        all_points.append(vector([x,y]))

    return node_count, cluster_count, edges_count, points, all_points

def read_solution_file(solution_file):
    edges_count = int(solution_file.readline())
    edges = []
    for i in range(edges_count):
        (start,end,weight) = solution_file.readline().split(' ')
        start = int(start)
        end = int(end)
        weight = float(weight)
        edges.append([start,end,weight])
    return edges_count, edges

if(len(sys.argv) < 4):
    print 'Usage: %s PROBLEM SOLUTION OUTPUT' % sys.argv[0]
    print 'Draw the nodes and clusters from file PROBLEM and the tour described by SOLUTION to OUTPUT.'
    print 'SOLUTION may be fractional. Acceptable OUTPUT file formats include PDF, EPS, SVG and PNG.'

else:
    problem_file = open(sys.argv[1], "r")
    solution_file = open(sys.argv[2], "r")

    (node_count, cluster_count, edges_count, points, all_points) = read_problem_file(problem_file)
    (edges_count, edges) = read_solution_file(solution_file)

    plot = list_plot([])

    max_x = max([p[0] for p in all_points])
    text_offset = vector([0,-1]) * max_x * 0.02

    print ('Drawing tour...')
    for k in range(edges_count):
        if edges[k][2] > 0.99:
            c = blue
        else:
            c = white.blend(red, 0.1 + 0.9 * edges[k][2])
        plot = plot + line([all_points[edges[k][0]], all_points[edges[k][1]]], color=c)

    print ('Drawing labels...')
    for i in range(node_count):
        plot = plot + text(str(i), all_points[i] + text_offset, color='gray')

    print ('Drawing clusters...')
    for i in range(cluster_count):
        plot = plot + list_plot(points[i], color='gray', figsize=FIGURE_SIZE,
                pointsize=POINT_SIZE)

        if(len(points[i]) > 1):
            vertices = chull_2d(list(set(points[i])))
            plot = plot + sum([line([vertices[k],vertices[(k+1)%len(vertices)]],
                    thickness=POINT_SIZE/20, color='gray') for k in range(len(vertices))])

    print("Writing file %s..." % sys.argv[3])
    save(plot, sys.argv[3])
