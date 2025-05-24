# Imports
import sys, math, random, matplotlib.pyplot as plt, networkx as nx, numpy as np
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QSizePolicy, QMessageBox
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.patches import Polygon
from scipy.spatial import Voronoi
from z3 import Solver, Bool, Or, Not, sat

# Main class
class GraphEditor(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ToC - Graph Builder")

        # Initializing the graph
        self.graph = nx.Graph()

        # Auxiliary variables
        self.k = 5

        self.node_size = 350
        self.node_radius = math.sqrt(self.node_size) * 2
        self.selected_node = None
        self.hovered_node = None
        self.node_colors = ['lightgray'] * len(self.graph.nodes)

        self.layout_timer = None
        self.animation_frame = 0
        self.animation_frames = 30
        self.interpolated_positions = {}

        self.mouseX = None
        self.mouseY = None

        # Layouts
        main_layout = QVBoxLayout()
        controls_layout = QHBoxLayout()

        # Expansion policy for the widgets
        widgets_exp_policy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        # Label for input form
        self.k_label = QLabel("Colors (k):")
        self.k_label.setSizePolicy(widgets_exp_policy)
        self.k_label.setMaximumWidth(250)

        # Input form for K
        self.k_input = QLineEdit()
        self.k_input.setSizePolicy(widgets_exp_policy)
        self.k_input.setMaximumWidth(250)
        self.k_input.setPlaceholderText(f"Enter k (Default = {self.k})")

        # Input button to solve the graph
        self.solve_button = QPushButton("Solve")
        self.solve_button.setSizePolicy(widgets_exp_policy)
        self.solve_button.setMaximumWidth(200)
        self.solve_button.clicked.connect(self.solve)

        # Input button to reset the graph
        self.reset_button = QPushButton("Reset")
        self.reset_button.setSizePolicy(widgets_exp_policy)
        self.reset_button.setMaximumWidth(200)
        self.reset_button.clicked.connect(self.reset_graph)

        # Input button to reset the graph
        self.spring_button = QPushButton("Reorganize")
        self.spring_button.setSizePolicy(widgets_exp_policy)
        self.spring_button.setMaximumWidth(200)
        self.spring_button.clicked.connect(self.animate_spring_layout)

        # Input button to generate the map
        self.map_button = QPushButton("Map")
        self.map_button.setSizePolicy(widgets_exp_policy)
        self.map_button.setMaximumWidth(200)
        self.map_button.clicked.connect(self.draw_voronoi_map)

        # Empty space taker
        dummy = QLabel("")
        dummy.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        controls_layout.addWidget(self.k_label)
        controls_layout.addWidget(self.k_input)
        controls_layout.addWidget(self.solve_button)
        controls_layout.addWidget(self.reset_button)
        controls_layout.addWidget(self.spring_button)
        # controls_layout.addWidget(self.map_button)
        controls_layout.addWidget(dummy)

        # Graph canvas
        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.figure)

        # Disabling autoscale
        self.ax.axis('off')

        # Handlers
        self.canvas.mpl_connect("button_press_event", self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)

        # Combine layouts
        main_layout.addLayout(controls_layout)
        main_layout.addWidget(self.canvas)
        self.setLayout(main_layout)

        self.draw_graph()

    # Handler of the on click event
    def on_click(self, event):
        # Verifying that the click is within the boundaries of the graph
        if not event.inaxes:
            return

        # Collecting the event data (i.e. position)
        x, y = event.xdata, event.ydata

        # Linear search to find the closest node
        clicked_node = self.get_node_at_position(x, y)

        # Left click
        if event.button == 1 and clicked_node is not None:
            # If this is the first selected node, cache it
            if self.selected_node is None:
                self.selected_node = clicked_node
            # Otherwise, manipulate the edges
            else :
                if self.graph.has_edge(self.selected_node, clicked_node):
                    self.graph.remove_edge(self.selected_node, clicked_node)
                else:
                    self.graph.add_edge(self.selected_node, clicked_node)
                self.selected_node = None
                self.reset_colors()

        # Right click
        if event.button == 3:
            # If there is a selected node, leave the selection
            if self.selected_node is not None:
                self.selected_node = None
            # Otherwise add/remove node
            else:
                if clicked_node is None:
                    new_node = len(self.graph)
                    self.graph.add_node(new_node, pos = (x, y))
                else:
                    self.graph.remove_node(clicked_node)
                    self.relabel_graph_sequentially()
                self.reset_colors()

        # Redrawing the graph
        self.draw_graph()

    # Auxiliary function that enlarges hovered nodes
    def on_hover(self, event):
        # Verifying that the click is within the boundaries of the graph
        if not event.inaxes:
            return

        # Updating the cached mouse position
        self.mouseX, self.mouseY = event.xdata, event.ydata

        # Looking for a node nearby
        self.hovered_node = self.get_node_at_position(self.mouseX, self.mouseY)

        # Redrawing the graph
        self.draw_graph()

    # Auxiliary function that find if the position overlaps with a node
    def get_node_at_position(self, x, y):
        # Initializing closest node
        closest = None
        shortest_sqr_distance = float('inf')

        # Transforming the given position in screen space
        screen_x, screen_y = self.ax.transData.transform((x, y))

        # Iterating all nodes, trying to find the closest one
        for node, (node_x, node_y) in nx.get_node_attributes(self.graph, 'pos').items():
            # Transforming the node position in screen space
            nscreen_x, nscreen_y = self.ax.transData.transform((node_x, node_y))

            # Computing the distance from node position and given position
            sqr_distance = (screen_x - nscreen_x) ** 2 + (screen_y - nscreen_y) ** 2

            # Verifying if the node is the closest
            if sqr_distance < self.node_radius ** 2 and sqr_distance < shortest_sqr_distance:
                closest = node
                shortest_sqr_distance = sqr_distance

        return closest

    # Auxiliary function that corrects the graph labels after an edge is deleted
    def relabel_graph_sequentially(self):
        # Extracting the old nodes
        old_nodes = sorted(self.graph.nodes())

        # Generating the mapping between the sets
        mapping = {}
        for new_label, old_label in enumerate(old_nodes):
            mapping[old_label] = new_label

        # Applying relabeling
        self.graph = nx.relabel_nodes(self.graph, mapping)

    # Auxiliary function that resets the colors of the graph
    def reset_colors(self):
        self.node_colors = ['lightgray'] * len(self.graph.nodes)

    # Auxiliary function that draw the graph based on current vertices and edges
    def draw_graph(self):
        # Clearing the board
        self.ax.clear()

        # Collecting the node positions
        positions = nx.get_node_attributes(self.graph, 'pos')

        # Initializing colors and sizes
        colors = ['lightgray'] * len(self.graph.nodes)
        sizes = [self.node_size] * len(self.graph.nodes)

        # Enlarging hovered node
        if self.hovered_node is not None and self.hovered_node in self.graph.nodes:
            sizes[self.hovered_node] = (self.node_size * 2)

        # Drawing a line from the currently selected node (if any)
        if self.selected_node is not None and self.mouseX is not None and self.mouseY is not None:
            self.ax.plot(
                [positions[self.selected_node][0], self.mouseX],
                [positions[self.selected_node][1], self.mouseY],
                color = 'gray',
                linestyle = '--',
                linewidth = 1.5,
            )

        # Drawing the graph
        nx.draw (
            self.graph, positions, ax=self.ax, with_labels=True,
            edge_color = 'black',
            node_color = self.node_colors, node_size = sizes)

        # Drawing the canvas
        self.ax.set_position([0, 0, 1, 1])
        self.figure.subplots_adjust(left=0, right=1, top=1, bottom=0)
        self.canvas.draw()

    # Auxiliary function that resets the graph
    def reset_graph(self):
        self.graph.clear()
        self.selected_node = None
        self.hovered_node = None
        self.draw_graph()

    # Auxiliary function that solve
    def solve(self):
        # Trying to parse the user input in the k_input to an int, otherwise use default value
        try:
            self.k = int(self.k_input.text()) if self.k_input.text() else 3
        except ValueError:
            self.k = 3

        # Initializing the SAT solver using Z3
        solver = Solver()

        # Extracting data from the current graph and initializing the colors range
        nodes = list(self.graph.nodes)
        colors = range(self.k)

        # Initializing the variables of the solver, in the form (vertex number, color number)
        variables = {(v, c): Bool(f"x({v},{c})") for v in nodes for c in colors}

        # Linear scan of the nodes, adding two clauses
        for n in nodes:
            # Clause: Every node must have a color
            solver.add(Or([variables[n, c] for c in colors]))

            # Clause: Every node has at most one color
            for c1 in colors:
                for c2 in colors:
                    if c1 < c2:
                        solver.add(Not(variables[n, c1]) | Not(variables[n, c2]))

        # Clause: Every neighbouring nodes must be colored differently
        for n1, n2 in self.graph.edges:
            for c in colors:
                solver.add(Not(variables[n1, c]) | Not(variables[n2, c]))

        # Verifying that the provided graph is SAT
        if solver.check() == sat:
            # Retrieving the solution
            model = solver.model()

            # Extrapolating each node's color from the solution and caching it
            for n in nodes:
                for c in colors:
                    if model.evaluate(variables[n, c]):
                        self.node_colors[n] = plt.cm.tab10(c % 10)
                        break

            # Drawing the graph
            self.draw_graph()
        else:
            QMessageBox.warning(self, "Unsolvable Graph", f"No valid coloring exists for k = {self.k}.")

    # Auxiliary function that animates the conversion of the graph to spring based layout by interpolating the positions
    def animate_spring_layout(self, duration_ms = 5000):
        # If the graph has only 2 nodes, nothing to do
        if len(self.graph.nodes) < 2:
            return

        # Getting current positions and computing target positions
        current_positions = nx.get_node_attributes(self.graph, 'pos')
        target_positions = nx.spring_layout(self.graph, seed=42)

        # Computing per-node delta
        for n in self.graph.nodes:
            x0, y0 = current_positions[n]
            x1, y1 = target_positions[n]
            self.interpolated_positions[n] = (x0, y0, x1, y1)

        # Computing the interval between each animation frame
        interval = duration_ms // self.animation_frames

        # Initializing and starting the animation
        self.layout_timer = QTimer(self)
        self.layout_timer.timeout.connect(self.step_layout_animation)
        self.layout_timer.start(interval)

    # Auxiliary function that execute an animation step for the spring layout
    def step_layout_animation(self):
        # Computing the interpolation value
        interpolation = self.animation_frame / self.animation_frames

        # Initializing the new positions
        new_positions = {}

        # Interpolating each node position
        for n, (x0, y0, x1, y1) in self.interpolated_positions.items():
            xt = x0 + (x1 - x0) * interpolation
            yt = y0 + (y1 - y0) * interpolation
            new_positions[n] = (xt, yt)

        # Setting the new positions and drawing the graph
        nx.set_node_attributes(self.graph, new_positions, 'pos')
        self.draw_graph()

        # Increasing the current frame
        self.animation_frame += 1

        # Verifying if the animation is complete
        if self.animation_frame > self.animation_frames:
            self.layout_timer.stop()

    # Auxiliary function that converts the graph into a map
    def draw_voronoi_map(self):
        # Clearing the space
        self.ax.clear()

        # Extracting the positions
        positions = nx.get_node_attributes(self.graph, 'pos')

        # Converting node positions to array
        node_indices = list(positions.keys())
        points = np.array([positions[n] for n in node_indices])
        vor = Voronoi(points)

        # Assign default or random color per node
        color_map = {n: random.randint(0, 9) for n in node_indices}

        # Draw Voronoi regions
        for i, region_index in enumerate(vor.point_region):
            region = vor.regions[region_index]
            if not -1 in region and len(region) > 0:
                polygon = [vor.vertices[j] for j in region]
                node = node_indices[i]
                color_id = color_map.get(node, 0)
                poly = Polygon(polygon, facecolor=plt.cm.tab10(color_id % 10), edgecolor='black', alpha=0.8)
                self.ax.add_patch(poly)

        self.ax.set_aspect('equal')
        self.ax.set_xlim(points[:, 0].min() - 0.1, points[:, 0].max() + 0.1)
        self.ax.set_ylim(points[:, 1].min() - 0.1, points[:, 1].max() + 0.1)
        self.ax.axis('off')
        self.canvas.draw()

# Main
if __name__ == '__main__':
    app = QApplication(sys.argv)
    editor = GraphEditor()
    editor.resize(1920, 1080)
    editor.show()
    sys.exit(app.exec_())