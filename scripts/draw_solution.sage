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

# data file
(node_count, cluster_count) = [ int(x) for x in raw_input().split(' ') ]

all_points = []
points = []
for i in range(cluster_count):
    points.append([])

for i in range(node_count):
    (x,y,cluster) = [ int(float(x)) for x in raw_input().split(' ') ]
    points[cluster].append((x,y))
    all_points.append((x,y))

# solutions file
edges_count = int(raw_input())
edges = []
for i in range(edges_count):
    edges.append([int(x) for x in raw_input().split(' ') ])


plot = list_plot([], xmax=100, xmin=0, ymax=100, ymin=0)

for i in range(cluster_count):
    plot = plot + list_plot(points[i], color='gray', figsize=FIGURE_SIZE,
            pointsize=POINT_SIZE)

    if(len(points[i]) > 1):
        vertices = chull_2d(list(set(points[i])))
        plot = plot + sum([line([vertices[k],vertices[(k+1)%len(vertices)]],
                thickness=POINT_SIZE/20, color='gray') for k in range(len(vertices))])

plot = plot + sum([line([all_points[edges[k][0]], all_points[edges[k][1]]]) for k in range(edges_count)])

save(plot, "tmp/gtsp.png")
