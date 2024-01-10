# flask --app website --debug run

import os
import subprocess
from flask import Flask, render_template, redirect, url_for, request, flash
from werkzeug.utils import secure_filename
from flask import session
from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, FileField, SubmitField, DecimalField, BooleanField, RadioField
from wtforms import Form, DecimalField, validators
from wtforms.validators import InputRequired, Optional, NumberRange
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('agg')
import numpy as np

app = Flask(__name__)
app.config['SECRET_KEY'] = 'Thisisasecretform'

# Input
UPLOAD_FOLDER = 'files'
ALLOWED_EXTENSIONS = {'stl', 'city.json', 'obj', 'ply'}
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
class UploadForm(FlaskForm): #Defines input parameters presented in User Interface I: Input file
    file = FileField("File", validators=[InputRequired()])
    street_distribution = BooleanField('Street distribution', validators=[Optional()], default=False)
    submit = SubmitField("Upload File")
class LoginForm(FlaskForm): #Defines input parameters presented in User Interface II: Input parameters

    # CFD
    d_flow = DecimalField('Flow direction:', validators=[Optional(), NumberRange(min=0, max=None)], default=90)
    h_user = DecimalField('Target evaluation height:', validators=[Optional(), NumberRange(min=0, max=None)], default=1.50)
    unit = DecimalField('Unit parameter:', validators=[Optional(), NumberRange(min=0, max=None)], default=1)
    BR = RadioField('Refinement boxes:', choices=[('no BR', 'no BR'), ('BR', 'BR'), ('BRL and BRH', 'BRL and BRH')], default='no BR')
    dim = RadioField('Dimension', choices=[('height', 'height'), ('dimension', 'dimension')], default='height')
    n_boxes = RadioField('Refinement boxes', choices=[('2', '2'), ('3', '3')], default='2')
    z0 = DecimalField('Maximum roughness length:', validators=[InputRequired(), NumberRange(min=0, max=None)])
    max = DecimalField('Maximum number of cells:', validators=[Optional(), NumberRange(min=0, max=None)], default=30e6, places=0)
    r_min = DecimalField('Maximum cell ratio:', validators=[Optional(), NumberRange(min=0, max=None)], default=1.20)
    d_street = DecimalField('Maximum building separation:', validators=[Optional()], default=2)
    target = BooleanField('Target building:', validators=[Optional()], default=False)
    distribution = BooleanField('Street distribution:', validators=[Optional()], default=False)
    x = DecimalField('x', validators=[Optional()])
    y = DecimalField('y', validators=[Optional()])

    f1_in = DecimalField('f1_in', validators=[Optional(), NumberRange(min=0, max=5)], default=0.5)
    f1_out = DecimalField('f1_out', validators=[Optional(), NumberRange(min=0, max=15)], default=0.5)
    f1_lat = DecimalField('f1_lat', validators=[Optional(), NumberRange(min=0, max=5)], default=0.5)
    f1_top = DecimalField('f1_top', validators=[Optional(), NumberRange(min=0, max=5)], default=1.5)

    f2_in = DecimalField('f2_in', validators=[Optional(), NumberRange(min=0, max=5)], default=2)
    f2_out = DecimalField('f2_out', validators=[Optional(), NumberRange(min=0, max=15)], default=6)
    f2_lat = DecimalField('f2_lat', validators=[Optional(), NumberRange(min=0, max=5)], default=2)
    f2_top = DecimalField('f2_top', validators=[Optional(), NumberRange(min=0, max=5)], default=2)

    f3_in = DecimalField('f3_in', validators=[Optional(), NumberRange(min=0, max=5)], default=3)
    f3_out = DecimalField('f3_out', validators=[Optional(), NumberRange(min=0, max=15)], default=9)
    f3_lat = DecimalField('f3_lat', validators=[Optional(), NumberRange(min=0, max=5)], default=3)
    f3_top = DecimalField('f3_top', validators=[Optional(), NumberRange(min=0, max=5)], default=3)

    # Geometric validations
    snaptolerance = DecimalField('Snap tolerance:', validators=[Optional(), NumberRange(min=0, max=None)], default=0.001, places=3)
    planaritytolerance = DecimalField('Planarity tolerance:', validators=[Optional(), NumberRange(min=0, max=None)], default=0.01)
    overlaptolerance = DecimalField('Overlap tolerance:', validators=[Optional(), NumberRange(min=-1, max=None)], default=-1)
    th_topo = DecimalField('Topological relationships threshold:', validators=[Optional()], default=0.5)
    th_slivertriangles = DecimalField('Sliver triangle threshold:', validators=[Optional(), NumberRange(min=0, max=None)], default=0.005)
    th_shortedges = DecimalField('Short edge threshold:', validators=[Optional(), NumberRange(min=0, max=None)], default=0.1)
    th_sharpangles = DecimalField('Sharp angle threshold:', validators=[Optional()], default=2)
    g_level = BooleanField('Ground level', validators=[Optional()], default=False)
    z_ground = DecimalField('Height:', validators=[Optional()])
    gs = BooleanField('Ground surfaces:', validators=[Optional()], default=False)
    h_gs = DecimalField('Maximum height:', validators=[Optional()], default=3)
    theta_gs = DecimalField('Maximum angle:', validators=[Optional()], default=45)

    submit = SubmitField("Upload File")

@app.route('/', methods=['GET', 'POST'])
def upload_file(): #Retrieves input from User Interface I: Input file
    form_file = UploadForm()
    global filename

    if form_file.validate_on_submit():
        file = request.files['file']
        session['street_distribution'] = form_file.street_distribution.data
        filename = request.files['file'].filename
        file.save(os.path.join(os.path.abspath(os.path.dirname(__file__)), app.config['UPLOAD_FOLDER'], secure_filename(file.filename)))
        return redirect("/form")

    return render_template("input.html", form=form_file)

def histogram(input_distances): #Creates histogram with street width distribution
    x = []
    with open(input_distances) as f:
        lines = f.readlines()
        for line in lines:
            if (float(line) <= 50):
                x.append(float(line))
    n, bins, patches = plt.hist(x, bins=25, facecolor='#62A39F')
    plt.title('Building separation distribution')
    plt.xlabel('Separations (m)')
    plt.xticks(np.arange(0, 50, 2))
    plt.ylabel('Frequency')
    max = n[0]
    idx = 0
    for i in range(len(n)):
        if n[i] > max:
            max = n[i]
            idx = i
    patches[idx].set_facecolor('#FFC300')
    print(idx)
    plt.savefig("website/static/output/histogram.png")

@app.route('/form', methods=['GET', 'POST'])
def form_file(): #Retrieves input parameters from users inserted in User Interface II: Input parameters
    name = "website/files/{}".format(filename)

    street_distribution = session.get('street_distribution')
    print(street_distribution)
    if street_distribution == 1:
        subprocess.run(['./website/histogram/cmake-build-debug/histogram', name, "website/static/output/distances.txt", "website/static/output/buildings.obj"])
        histogram('website/static/output/distances.txt')

    form_parameters = LoginForm()
    if form_parameters.validate_on_submit():
        session['d_flow'] = form_parameters.d_flow.data
        session['h_user'] = form_parameters.h_user.data
        session['unit'] = form_parameters.unit.data
        session['BR'] = form_parameters.BR.data
        session['dim'] = form_parameters.dim.data
        session['n_boxes'] = form_parameters.n_boxes.data

        session['f1_in'] = form_parameters.f1_in.data
        session['f1_out'] = form_parameters.f1_out.data
        session['f1_lat'] = form_parameters.f1_lat.data
        session['f1_top'] = form_parameters.f1_top.data

        session['f2_in'] = form_parameters.f2_in.data
        session['f2_out'] = form_parameters.f2_out.data
        session['f2_lat'] = form_parameters.f2_lat.data
        session['f2_top'] = form_parameters.f2_top.data

        session['f3_in'] = form_parameters.f3_in.data
        session['f3_out'] = form_parameters.f3_out.data
        session['f3_lat'] = form_parameters.f3_lat.data
        session['f3_top'] = form_parameters.f3_top.data

        session['z0'] = form_parameters.z0.data
        session['max'] = form_parameters.max.data
        session['r_min'] = form_parameters.r_min.data
        session['d_street'] = form_parameters.d_street.data
        session['distribution'] = form_parameters.distribution.data
        session['target'] = form_parameters.target.data
        session['x'] = form_parameters.x.data
        session['y'] = form_parameters.y.data

        session['snaptolerance'] = form_parameters.snaptolerance.data
        session['planaritytolerance'] = form_parameters.planaritytolerance.data
        session['overlaptolerance'] = form_parameters.overlaptolerance.data

        session['th_topo'] = form_parameters.th_topo.data
        session['g_level'] = form_parameters.g_level.data
        session['z_ground'] = form_parameters.z_ground.data
        session['th_slivertriangles'] = form_parameters.th_slivertriangles.data
        session['th_shortedges'] = form_parameters.th_shortedges.data
        session['th_sharpangles'] = form_parameters.th_sharpangles.data

        session['gs'] = form_parameters.gs.data
        session['h_gs'] = form_parameters.h_gs.data
        session['theta_gs'] = form_parameters.theta_gs.data

        return redirect("/home")

    return render_template("form.html", form=form_parameters, street_distribution=street_distribution)
def user_parameters(): #Creates a file with user-defined parameters used to run the C++ code "backend"
    wd = session.get('d_flow')
    n_boxes = session.get('n_boxes')
    unit = session.get('unit')
    h_user = session.get('h_user')
    dim = session.get('dim')
    z0 = session.get('z0')
    N_max = session.get('max')
    r_min = session.get('r_min')
    d_separation = session.get('d_street')
    BR = session.get('BR')
    target = session.get('target')
    x = session.get('x')
    y = session.get('y')

    f1_in = session.get('f1_in')
    f1_out = session.get('f1_out')
    f1_lat = session.get('f1_lat')
    f1_top = session.get('f1_top')

    f2_in = session.get('f2_in')
    f2_out = session.get('f2_out')
    f2_lat = session.get('f2_lat')
    f2_top = session.get('f2_top')

    f3_in = session.get('f3_in')
    f3_out = session.get('f3_out')
    f3_lat = session.get('f3_lat')
    f3_top = session.get('f3_top')

    th_topo = session.get('th_topo')
    ground_level = session.get('g_level')
    z_ground = session.get('z_ground')

    th_slivers = session.get('th_slivertriangles')
    th_angles = session.get('th_sharpangles')
    th_edges = session.get('th_shortedges')

    gs = session.get('gs')
    h_gs = session.get('h_gs')
    theta_gs = session.get('theta_gs')

    snaptolerance = session.get('snaptolerance')
    planaritytolerance = session.get('planaritytolerance')
    overlaptolerance = session.get('overlaptolerance')

    if BR == 'no BR':
        BR_option = 0
    if BR == 'BR':
        BR_option = 1
    if BR == 'BRL and BRH':
        BR_option = 2

    if dim == 'height':
        dim_option = 0
    if dim == 'dimension':
        dim_option = 1

    if n_boxes == '2':
        n_option = 2
    if n_boxes == '3':
        n_option = 3

    with open('website/files/user_parameters.txt', 'w') as f:
        f.write('name			{0}\n'.format(filename))
        f.write('buildings              website/static/output/buildings.obj\n')
        f.write('wd                     {0}\n'.format(wd))
        f.write('unit                   {0}\n'.format(unit))
        f.write('n_boxes                {0}\n'.format(n_boxes))
        f.write('h_user                 {0}\n'.format(h_user))
        f.write('z0                     {0}\n'.format(z0))
        f.write('N_max                  {0}\n'.format(N_max))
        f.write('r_min                  {0}\n'.format(r_min))
        f.write('d_separation           {0}\n'.format(d_separation))
        f.write('BR                     {0}\n'.format(BR_option))
        f.write('dim                    {0}\n'.format(dim_option))
        f.write('target                 {0}\n'.format(int(target)))
        f.write('x                      {0}\n'.format(x))
        f.write('y                      {0}\n'.format(y))
        f.write('f1                     {} {} {} {}\n'.format(f1_in, f1_out, f1_lat, f1_top))
        f.write('f2                     {} {} {} {}\n'.format(f2_in, f2_out, f2_lat, f2_top))
        f.write('f3                     {} {} {} {}\n'.format(f3_in, f3_out, f3_lat, f3_top))
        f.write('snap_tolerance         {0}\n'.format(snaptolerance))
        f.write('planarity_tolerance    {0}\n'.format(planaritytolerance))
        f.write('overlap_tolerance      {0}\n'.format(overlaptolerance))
        f.write('th_topo                {0}\n'.format(th_topo))
        f.write('ground_level           {0}\n'.format(int(ground_level)))
        f.write('z_ground               {0}\n'.format(z_ground))
        f.write('th_slivers             {0}\n'.format(th_slivers))
        f.write('th_angles              {0}\n'.format(th_angles))
        f.write('th_edges               {0}\n'.format(th_edges))

        if gs == True:
            f.write('h_gs                   {0}\n'.format(h_gs))
            f.write('theta_gs               {0}\n'.format(theta_gs))
        else:
            f.write('h_gs                   3\n')
            f.write('theta_gs               45\n')
        f.write('save_gs                website/static/output/gs.obj\n')
        f.write('topocheck              website/static/output/topocheck.obj                     website/static/output/topocheck.txt\n')
        f.write('slivers                website/static/output/sliver_triangles.obj              website/static/output/sliver_triangles.txt\n')
        f.write('sharpangles            website/static/output/sharp_angles.obj                  website/static/output/sharp_angles.txt\n')
        f.write('shortedges             website/static/output/short_edges.obj                   website/static/output/short_edges.txt\n')
        f.write('overlappingbuildings   website/static/output/overlapping_buildings.obj         website/static/output/overlapping_buildings.txt\n')
        f.write('val3dity               website/static/output/report_val3dity.json\n')
        f.write('rotated_model          website/static/output/rotated_model.obj\n')
        f.write('blockMeshDict          website/backend/data/blockMeshDict   website/static/output/blockMeshDict                     website/static/output/blockMeshDict_noNmax              website/static/output/blockMeshDict_Nmax\n')
        f.write('snappyHexMeshDict	website/backend/data/snappyHexMeshDict{0}      website/static/output/snappyHexMeshDict                 website/static/output/snappyHexMeshDict_noNmax          website/static/output/snappyHexMeshDict_Nmax\n'.format(n_option))
        f.write('domain                 website/static/output/domain.obj                        website/static/output/domain_noNmax.obj                 website/static/output/domain_Nmax.obj\n')
        f.write('b_volume               website/static/output/invalid_building_volumes.obj	website/static/output/invalid_building_volumes_noNmax.obj	website/static/output/invalid_building_volumes_Nmax.obj\n')
        f.write('b_separations          website/static/output/b_separations_buildings.obj       website/static/output/b_separations_edges.obj		website/static/output/b_separations_buildings_noNmax.obj       website/static/output/b_separations_edges_noNmax.obj	website/static/output/b_separations_buildings_Nmax.obj       website/static/output/b_separations_edges_Nmax.obj\n')
        f.write('RoI                    website/static/output/RoI.obj                           website/static/output/RoI.txt                           website/static/output/RoI_cylinder.obj\n')
        f.write('output_for_UI   	website/static/output/output_for_UI.txt\n')
def read_output_UI(): #Reads results from C++ code "backend"
    with open('website/static/output/output_for_UI.txt') as f:
        lines = f.readlines()

        for line in lines:
            l = line.split()
            if l[0] == "valid":
                global valid
                valid = int(l[1])
            if l[0] == "n_topo":
                global n_topo
                n_topo = int(l[1])
            if l[0] == "n_sharp_angles":
                global n_sharp_angles
                n_sharp_angles = int(l[1])
            if l[0] == "n_short_edges":
                global n_short_edges
                n_short_edges = int(l[1])
            if l[0] == "n_slivers":
                global n_slivers
                n_slivers = int(l[1])
            if l[0] == "cell":
                global h_cell
                global w_cell
                h_cell = float(l[1])
                w_cell = float(l[2])
            if l[0] == "cell_GR_noNmax":
                global h_cell_GR_noNmax
                global w_cell_GR_noNmax
                h_cell_GR_noNmax = float(l[1])
                w_cell_GR_noNmax = float(l[2])
            if l[0] == "cell_GR_Nmax":
                global h_cell_GR_Nmax
                global w_cell_GR_Nmax
                h_cell_GR_Nmax = float(l[1])
                w_cell_GR_Nmax = float(l[2])
            if l[0] == "eh":
                global eh1a
                eh1a = float(l[1])
                global eh1b
                eh1b = float(l[2])
            if l[0] == "eh_GR_noNmax":
                global eh_GR_noNmax
                eh_GR_noNmax = float(l[1])
            if l[0] == "eh_GR_Nmax":
                global eh_GR_Nmax
                eh_GR_Nmax = float(l[1])
            if l[0] == "n_RoI":
                global n_RoI
                n_RoI = int(l[1])
            if l[0] == "per_RoI":
                global per_RoI
                per_RoI = float(l[1])
            if l[0] == "BV":
                global BV
                BV = float(l[1])
            if l[0] == "BV_valid":
                global BV_valid
                BV_valid = float(l[1])
            if l[0] == "BV_valid_GR_noNmax":
                global BV_valid_GR_noNmax
                BV_valid_GR_noNmax = float(l[1])
            if l[0] == "BV_valid_GR_Nmax":
                global BV_valid_GR_Nmax
                BV_valid_GR_Nmax = float(l[1])
            if l[0] == "S_valid":
                global S_valid
                S_valid = float(l[1])
            if l[0] == "S_valid_GR_noNmax":
                global S_valid_GR_noNmax
                S_valid_GR_noNmax = float(l[1])
            if l[0] == "S_valid_GR_Nmax":
                global S_valid_GR_Nmax
                S_valid_GR_Nmax = float(l[1])
            if l[0] == "N":
                global N
                N = float(l[1])
            if l[0] == "N_GR_noNmax":
                global N_GR_noNmax
                N_GR_noNmax = float(l[1])
            if l[0] == "N_GR_Nmax":
                global N_GR_Nmax
                N_GR_Nmax = float(l[1])
            if l[0] == "BR":
                global BR
                BR = float(l[1])
            if l[0] == "BRL":
                global BRL
                BRL = float(l[1])
            if l[0] == "BRH":
                global BRH
                BRH = float(l[1])
            if l[0] == "hx":
                global hx
                hx = float(l[1])
            if l[0] == "hy":
                global hy
                hy = float(l[1])
            if l[0] == "hz":
                global hz
                hz = float(l[1])
            if l[0] == "BR_noNmax":
                global BR_noNmax
                BR_noNmax = float(l[1])
            if l[0] == "BRL_noNmax":
                global BRL_noNmax
                BRL_noNmax = float(l[1])
            if l[0] == "BRH_noNmax":
                global BRH_noNmax
                BRH_noNmax = float(l[1])
            if l[0] == "hx_noNmax":
                global hx_noNmax
                hx_noNmax = float(l[1])
            if l[0] == "hy_noNmax":
                global hy_noNmax
                hy_noNmax = float(l[1])
            if l[0] == "hz_noNmax":
                global hz_noNmax
                hz_noNmax = float(l[1])
            if l[0] == "BR_Nmax":
                global BR_Nmax
                BR_Nmax = float(l[1])
            if l[0] == "BRL_Nmax":
                global BRL_Nmax
                BRL_Nmax = float(l[1])
            if l[0] == "BRH_Nmax":
                global BRH_Nmax
                BRH_Nmax = float(l[1])
            if l[0] == "hx_noNmax":
                global hx_Nmax
                hx_Nmax = float(l[1])
            if l[0] == "hy_Nmax":
                global hy_Nmax
                hy_Nmax = float(l[1])
            if l[0] == "hz_Nmax":
                global hz_Nmax
                hz_Nmax = float(l[1])
            if l[0] == "overlap":
                global n_overlap
                n_overlap = int(l[1])
            if l[0] == "n_topo_errors":
                global n_topo_errors
                n_topo_errors = int(l[1])
            if l[0] == "n_topo_warnings":
                global n_topo_warnings
                n_topo_warnings = int(l[1])
            if l[0] == "dmax":
                global dmax
                dmax = float(l[1])

@app.route("/home", methods=['GET'])
def table(): #Provides results to users in User Interface III: Results
    user_parameters()
    name = "website/files/{}".format(filename)
    subprocess.run(['./website/backend/cmake-build-debug/backend',  name, "website/files/user_parameters.txt"])
    read_output_UI()
    headings = ("Content", " ", "Details", "Output files")

    topo = n_topo
    topo_errors = n_topo_errors
    topo_warnings = n_topo_warnings
    sharp_angles = n_sharp_angles
    short_edges = n_short_edges
    slivers = n_slivers
    overlap = n_overlap
    hcell = round(h_cell, 2)
    wcell = round(w_cell, 2)
    hcell2 = round(h_cell_GR_noNmax, 2)
    wcell2 = round(w_cell_GR_noNmax, 2)
    hcell3 = round(h_cell_GR_Nmax, 2)
    wcell3 = round(w_cell_GR_Nmax, 2)
    eh_a = round(eh1a, 2)
    eh_b = round(eh1b, 2)
    eh2 = round(eh_GR_noNmax, 2)
    eh3 = round(eh_GR_Nmax, 2)
    RoI1 = n_RoI
    RoI2 = round(per_RoI, 2)
    n_BV_valid = round(BV_valid, 2)
    n_BV_valid1 = round(BV_valid_GR_noNmax, 2)
    n_BV_valid2 = round(BV_valid_GR_Nmax, 2)
    n_S_valid = round(S_valid, 2)
    n_S_valid1 = round(S_valid_GR_noNmax, 2)
    n_S_valid2 = round(S_valid_GR_Nmax, 2)
    direction = session['d_flow']
    N1 = round(N)
    N2 = round(N_GR_noNmax)
    N3 = round(N_GR_Nmax)
    BR1 = round(BR, 2)
    BRH1 = round(BRH, 2)
    BRL1 = round(BRL, 2)
    hx1 = round(hx, 2)
    hy1 = round(hy, 2)
    hz1 = round(hz, 2)
    BR2 = round(BR_noNmax, 2)
    BRH2 = round(BRH_noNmax, 2)
    BRL2 = round(BRL_noNmax, 2)
    hx2 = round(hx_noNmax, 2)
    hy2 = round(hy_noNmax, 2)
    hz2 = round(hz_noNmax, 2)
    BR3 = round(BR_Nmax, 2)
    BRH3 = round(BRH_Nmax, 2)
    BRL3 = round(BRL_Nmax, 2)
    hx3 = round(hx_Nmax, 2)
    hy3 = round(hy_Nmax, 2)
    hz3 = round(hz_Nmax, 2)
    hmax = round(dmax, 2)

    if valid == 0:
        check_validity = 'INVALID :('
    else:
        check_validity = 'VALID :)'

    return render_template("output.html", headings=headings, validity=check_validity, topo=n_topo,topo_errors=topo_errors, topo_warnings = topo_warnings,
                           sharp_angles=sharp_angles, short_edges=short_edges, slivers=slivers, hcell=hcell,
                           wcell=wcell, hcell2=hcell2, wcell2=wcell2, hcell3=hcell3, wcell3=wcell3, eh_a=eh_a,
                           eh_b=eh_b, eh2=eh2, eh3=eh3, RoI1=RoI1, RoI2=RoI2, n_BV_valid=n_BV_valid,
                           n_BV_valid1=n_BV_valid1, n_BV_valid2=n_BV_valid2, n_S_valid=n_S_valid, n_S_valid1=n_S_valid1,
                           n_S_valid2=n_S_valid2, direction=direction, N1=N1, N2=N2, N3=N3, BR1=BR1, BRH1=BRH1,
                           BRL1=BRL1, hx1=hx1, hy1=hy1, hz1=hz1, BR2=BR2, BRH2=BRH2, BRL2=BRL2, hx2=hx2, hy2=hy2,
                           hz2=hz2, BR3=BR3, BRH3=BRH3, BRL3=BRL3, hx3=hx3, hy3=hy3, hz3=hz3, overlap=overlap, h_max=hmax)

@app.route("/explanation", methods=['GET'])
def explanation(): #Renders User Interface IV: Explanations
    return render_template("explanation.html")

if __name__ == '__main__':
    user_parameters()
