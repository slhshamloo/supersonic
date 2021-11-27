import math
import numpy as np
from matplotlib import pyplot as plt


def double_solve(func1, func2, tn, y10=0, y20=0, t0=0, step=0.01):
    n = int((tn - t0) / step)
    t = np.linspace(t0, tn, n)
    h = step
    if len(t) > 1:
        h = t[1] - t[0]
    y1 = np.empty(n)
    y2 = np.empty(n)
    y1[0] = y10
    y2[0] = y20

    for i in range(1, n):
        k11 = func1(t[i-1], y1[i-1], y2[i-1])
        k21 = func2(t[i-1], y1[i-1], y2[i-1])
        k12 = func1(t[i-1] + h/2, y1[i-1] + h/2 * k11, y2[i-1] + h/2 * k21)
        k22 = func2(t[i-1] + h/2, y1[i-1] + h/2 * k11, y2[i-1] + h/2 * k21)
        k13 = func1(t[i-1] + h/2, y1[i-1] + h/2 * k12, y2[i-1] + h/2 * k22)
        k23 = func2(t[i-1] + h/2, y1[i-1] + h/2 * k12, y2[i-1] + h/2 * k22)
        k14 = func1(t[i], y1[i-1] + h * k13, y2[i-1] + h * k23)
        k24 = func2(t[i], y1[i-1] + h * k13, y2[i-1] + h * k23)
        y1[i] = y1[i-1] + (k11 + 2*k12 + 2*k13 + k14) * h / 6
        y2[i] = y2[i-1] + (k21 + 2*k22 + 2*k23 + k24) * h / 6
    
    return t, y1, y2


def integrate(der, tn, y0=0, t0=0, step=0.01):
    n = int((tn - t0) / step)
    t = np.linspace(t0, tn, n)
    h = step
    if len(t) > 1:
        h = t[1] - t[0]
    y = np.empty(n)
    y[0] = y0

    for i in range(1, n):
        y[i] = y[i-1] + h*der[i-1]
    return t, y


def cw(v):
    if v < 240:
        return 0.2
    if v < 275:
        return 0.2 + (v-240)/3500
    if v < 377:
        return 0.21 + (v-275)/102 * 0.165
    if v < 439:
        return 0.375 + (v-377)/62 * 0.025
    if v < 539:
        return 0.4
    if v < 686:
        return 0.4 - (v-539)/147 * 0.05
    else:
        return 0.35 - (v-686)/343 * 0.08


def ax(t, vx, vy):
    return - 0.001255*math.pi/6 * cw(math.sqrt(vx**2 + vy**2)) * math.sqrt(vx**2 + vy**2)*vx


def ay(t, vx, vy):
    return - 0.001255*math.pi/6 * cw(math.sqrt(vx**2 + vy**2)) * math.sqrt(vx**2 + vy**2)*vy - 9.8


def axr(t, vx, vy):
    return - 0.01255*math.pi/4 * cw(math.sqrt(vx**2 + vy**2)) * math.sqrt(vx**2 + vy**2)*vx


def ayr(t, vx, vy):
    return - 0.01255*math.pi/4 * cw(math.sqrt(vx**2 + vy**2)) * math.sqrt(vx**2 + vy**2)*vy - 9.8


def axl(t, vx, vy):
    return - 5*0.01255*math.pi * cw(math.sqrt(vx**2 + vy**2)) * math.sqrt(vx**2 + vy**2)*vx


def ayl(t, vx, vy):
    return - 5*0.01255*math.pi * cw(math.sqrt(vx**2 + vy**2)) * math.sqrt(vx**2 + vy**2)*vy - 9.8


def y_trim(t, x, y):
    tt = []
    xt = []
    yt = []
    for i in range(len(t)):
        if y[i] >= 0:
            tt.append(t[i])
            xt.append(x[i])
            yt.append(y[i])
        else:
            break
    return tt, xt, yt


def v_trim(t, vx, vy, x, y):
    tt = []
    vxt = []
    vyt = []
    xt = []
    yt = []
    for i in range(len(t)):
        if y[i] >= 0:
            tt.append(t[i])
            vxt.append(vx[i])
            vyt.append(vy[i])
            xt.append(x[i])
            yt.append(y[i])
        else:
            break
    return tt, vxt, vyt, xt, yt


def coeff():
    plt.figure(figsize=(12.0, 8.0))
    plt.plot([0, 240, 275, 377, 439, 539, 686, 1029],
             [0.2, 0.2, 0.21, 0.375, 0.4, 0.4, 0.35, 0.27])
    plt.xlabel('$v$ $(\mathrm{m}/\mathrm{s})$', fontsize=20)
    plt.ylabel('$C_w$', fontsize=20)
    plt.grid('True')
    plt.savefig('Cw.pdf', bbox_inches='tight')
    plt.show()

def xy():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(ax, ay, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)

    t, vxr, vyr = double_solve(axr, ayr, 100, v0, v0)
    t, xr = integrate(vxr, tn=100)
    t, yr = integrate(vyr, tn=100)

    x0 = v0 * t
    y0 = -4.9 * t**2 + v0 * t

    tt, xt, yt = y_trim(t, x, y)
    trt, xrt, yrt = y_trim(t, xr, yr)
    t0t, x0t, y0t = y_trim(t, x0, y0)

    plt.figure(figsize=(12.0, 8.0))
    plt.axes().set_aspect('equal')
    plt.grid('True')

    plt.plot(x0t, y0t, 'r', label='without drag')
    plt.plot(xt, yt, 'g', label='with drag, 30kg ball')
    plt.plot(xrt, yrt, 'b', label='with drag, 2kg ball')

    plt.xlabel('$x (\mathrm{m})$', fontsize='15')
    plt.ylabel('$y (\mathrm{m})$', fontsize='15')
    plt.legend()

    plt.savefig('xy.pdf', bbox_inches='tight')
    plt.show()


def yplot():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(ax, ay, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)

    t, vxr, vyr = double_solve(axr, ayr, 100, v0, v0)
    t, xr = integrate(vxr, tn=100)
    t, yr = integrate(vyr, tn=100)

    x0 = v0 * t
    y0 = -4.9 * t**2 + v0 * t

    tt, xt, yt = y_trim(t, x, y)
    trt, xrt, yrt = y_trim(t, xr, yr)
    t0t, x0t, y0t = y_trim(t, x0, y0)

    plt.figure(figsize=(12.0, 8.0))
    plt.grid('True')

    plt.plot(t0t, y0t, 'r', label='without drag')
    plt.plot(tt, yt, 'g', label='with drag, 30kg ball')
    plt.plot(trt, yrt, 'b', label='with drag, 2kg ball')

    plt.xlabel('$t (\mathrm{s})$', fontsize='20')
    plt.ylabel('$y (\mathrm{m})$', fontsize='20')
    plt.legend()

    plt.savefig('yplot.pdf', bbox_inches='tight')
    plt.show()


def xplot():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(ax, ay, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)

    t, vxr, vyr = double_solve(axr, ayr, 100, v0, v0)
    t, xr = integrate(vxr, tn=100)
    t, yr = integrate(vyr, tn=100)

    x0 = v0 * t
    y0 = -4.9 * t**2 + v0 * t

    tt, xt, yt = y_trim(t, x, y)
    trt, xrt, yrt = y_trim(t, xr, yr)
    t0t, x0t, y0t = y_trim(t, x0, y0)

    plt.figure(figsize=(12.0, 8.0))
    plt.grid('True')

    plt.plot(t0t, x0t, 'r', label='without drag')
    plt.plot(tt, xt, 'g', label='with drag, 30kg ball')
    plt.plot(trt, xrt, 'b', label='with drag, 2kg ball')

    plt.xlabel('$t (\mathrm{s})$', fontsize='20')
    plt.ylabel('$x (\mathrm{m})$', fontsize='20')
    plt.legend()

    plt.savefig('xplot.pdf', bbox_inches='tight')
    plt.show()


def vyplot():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(ax, ay, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)

    t, vxr, vyr = double_solve(axr, ayr, 100, v0, v0)
    t, xr = integrate(vxr, tn=100)
    t, yr = integrate(vyr, tn=100)

    x0 = v0 * t
    y0 = -4.9 * t**2 + v0 * t
    vx0 = v0 * np.ones(10000)
    vy0 = -9.8 * t + v0

    tt, vxt, vyt, xt, yt = v_trim(t, vx, vy, x, y)
    trt, vxrt, vyrt, xrt, yrt = v_trim(t, vxr, vyr, xr, yr)
    t0t, vx0t, vy0t, x0t, y0t = v_trim(t, vx0, vy0, xr, y0)

    plt.figure(figsize=(12.0, 8.0))
    plt.grid('True')

    plt.plot(t0t, vy0t, 'r', label='without drag')
    plt.plot(tt, vyt, 'g', label='with drag, 30kg ball')
    plt.plot(trt, vyrt, 'b', label='with drag, 2kg ball')

    plt.xlabel('$t (\mathrm{s})$', fontsize='20')
    plt.ylabel('$\dot{y} (\mathrm{m}/\mathrm{s})$', fontsize='20')
    plt.legend()

    plt.savefig('vyplot.pdf', bbox_inches='tight')
    plt.show()


def vxplot():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(ax, ay, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)

    t, vxr, vyr = double_solve(axr, ayr, 100, v0, v0)
    t, xr = integrate(vxr, tn=100)
    t, yr = integrate(vyr, tn=100)

    x0 = v0 * t
    y0 = -4.9 * t**2 + v0 * t
    vx0 = v0 * np.ones(10000)
    vy0 = -9.8 * t + v0

    tt, vxt, vyt, xt, yt = v_trim(t, vx, vy, x, y)
    trt, vxrt, vyrt, xrt, yrt = v_trim(t, vxr, vyr, xr, yr)
    t0t, vx0t, vy0t, x0t, y0t = v_trim(t, vx0, vy0, xr, y0)

    plt.figure(figsize=(12.0, 8.0))
    plt.grid('True')

    plt.plot(t0t, vx0t, 'r', label='without drag')
    plt.plot(tt, vxt, 'g', label='with drag, 30kg ball')
    plt.plot(trt, vxrt, 'b', label='with drag, 2kg ball')

    plt.xlabel('$t (\mathrm{s})$', fontsize='20')
    plt.ylabel('$\dot{x} (\mathrm{m}/\mathrm{s})$', fontsize='20')
    plt.legend()

    plt.savefig('vxplot.pdf', bbox_inches='tight')
    plt.show()


def yvy():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(ax, ay, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)

    t, vxr, vyr = double_solve(axr, ayr, 100, v0, v0)
    t, xr = integrate(vxr, tn=100)
    t, yr = integrate(vyr, tn=100)

    x0 = v0 * t
    y0 = -4.9 * t**2 + v0 * t
    vx0 = v0 * np.ones(10000)
    vy0 = -9.8 * t + v0

    tt, vxt, vyt, xt, yt = v_trim(t, vx, vy, x, y)
    trt, vxrt, vyrt, xrt, yrt = v_trim(t, vxr, vyr, xr, yr)
    t0t, vx0t, vy0t, x0t, y0t = v_trim(t, vx0, vy0, x0, y0)

    plt.figure(figsize=(12.0, 8.0))
    plt.grid('True')

    plt.plot(y0t, vy0t, 'r', label='without drag')
    plt.plot(yt, vyt, 'g', label='with drag, 30kg ball')
    plt.plot(yrt, vyrt, 'b', label='with drag, 2kg ball')

    plt.xlabel('$y (\mathrm{m})$', fontsize='20')
    plt.ylabel('$\dot{y} (\mathrm{m}/\mathrm{s})$', fontsize='20')
    plt.legend()

    plt.savefig('yvy.pdf', bbox_inches='tight')
    plt.show()


def xvx():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(ax, ay, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)

    t, vxr, vyr = double_solve(axr, ayr, 100, v0, v0)
    t, xr = integrate(vxr, tn=100)
    t, yr = integrate(vyr, tn=100)

    x0 = v0 * t
    y0 = -4.9 * t**2 + v0 * t
    vx0 = v0 * np.ones(10000)
    vy0 = -9.8 * t + v0

    tt, vxt, vyt, xt, yt = v_trim(t, vx, vy, x, y)
    trt, vxrt, vyrt, xrt, yrt = v_trim(t, vxr, vyr, xr, yr)
    t0t, vx0t, vy0t, x0t, y0t = v_trim(t, vx0, vy0, x0, y0)

    plt.figure(figsize=(12.0, 8.0))
    plt.grid('True')

    plt.plot(x0t, vx0t, 'r', label='without drag')
    plt.plot(xt, vxt, 'g', label='with drag, 30kg ball')
    plt.plot(xrt, vxrt, 'b', label='with drag, 2kg ball')

    plt.xlabel('$x (\mathrm{m})$', fontsize='20')
    plt.ylabel('$\dot{x} (\mathrm{m}/\mathrm{s})$', fontsize='20')
    plt.legend()

    plt.savefig('xvx.pdf', bbox_inches='tight')
    plt.show()


def light():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(axl, ayl, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)
    tt, xt, yt = y_trim(t, x, y)

    plt.figure(figsize=(12.0, 8.0))
    plt.axes().set_aspect('equal')
    plt.grid('True')

    plt.plot(xt, yt)
    plt.vlines(x=x[np.argmax(y)], ymin = 0, ymax=np.max(y), linestyles='--', colors='r')

    plt.xlabel('$x (\mathrm{m})$', fontsize='15')
    plt.ylabel('$y (\mathrm{m})$', fontsize='15')
    plt.savefig('light.pdf', bbox_inches='tight')
    plt.show()


def heavy():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(ax, ay, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)
    tt, xt, yt = y_trim(t, x, y)

    plt.figure(figsize=(12.0, 8.0))
    plt.axes().set_aspect('equal')
    plt.grid('True')

    plt.plot(xt, yt)
    plt.vlines(x=x[np.argmax(y)], ymin = 0, ymax=np.max(y), linestyles='--', colors='r')

    plt.xlabel('$x (\mathrm{m})$', fontsize='15')
    plt.ylabel('$y (\mathrm{m})$', fontsize='15')
    plt.savefig('heavy.pdf', bbox_inches='tight')
    plt.show()


def num():
    v0 = 600 / math.sqrt(2)

    t, vx, vy = double_solve(ax, ay, 100, v0, v0)
    t, x = integrate(vx, tn=100)
    t, y = integrate(vy, tn=100)

    t, vxr, vyr = double_solve(axr, ayr, 100, v0, v0)
    t, xr = integrate(vxr, tn=100)
    t, yr = integrate(vyr, tn=100)

    x0 = v0 * t
    y0 = -4.9 * t**2 + v0 * t

    tt, xt, yt = y_trim(t, x, y)
    trt, xrt, yrt = y_trim(t, xr, yr)
    t0t, x0t, y0t = y_trim(t, x0, y0)

    print(np.max(tt), np.max(xt), np.max(yt))
    print(np.max(trt), np.max(xrt), np.max(yrt))
    print(np.max(t0t), np.max(x0t), np.max(y0t))


if __name__ == '__main__':
    yvy()
