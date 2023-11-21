from django.shortcuts import render
from sympy.calculus.util import continuous_domain
import pandas as pd
import numpy as np
from sympy import *
import math
import matplotlib.pyplot as plt
from sympy.plotting import plot
from scipy.interpolate import interp1d

def menu_page(request):
    return render(request,template_name='menu.html')

def biseccion_page(request):
    if request.method == 'POST':
        x, y, z = symbols('x y z')
        Fun = request.POST.get('funcion')
        a = float(request.POST.get('a'))
        b = float(request.POST.get('b'))
        tipo_error = int(request.POST.get('tipo_de_tolerancia'))
        num_tol = float(request.POST.get('numero_de_tolerancia'))
        Niter = int(request.POST.get('iteraciones'))
        inia = a
        inib = b
        
        df=0

        if a > b:
            temp = a
            a, b = b, temp

        if tipo_error == 1: 
            tipo = "decimales correctos"
            Tol = 0.5*(10**-num_tol)
        elif tipo_error == 2: 
            tipo = "cifras significativas"
            Tol = 5*(10**-num_tol)
        
        Fun = sympify(Fun)
        func = lambda m: Fun.evalf(15, subs = {x: m})
        dom = continuous_domain(Fun, x, S.Reals)
        interval = Interval(a,b)

        if not interval.is_subset(dom):
            mensaje = f"La función no es continua en el intervalo [{a},{b}]. El método falla."
            graph = plot(Fun, xlabel='x', ylabel='y', show=False)
            graph.save('Grafica.png')
            context={'mensaje':mensaje,'df':df,'fun':Fun,'a':inia,'b':inib,'tipo_error':tipo,'num_tol':num_tol,'niter':Niter,'grafica':'../Grafica.png'}
            return render(request,template_name='1-biseccion.html',context=context)

        x_vals, fm, E, iters = [], [], [], []
        fa = func(a)
        fb = func(b)
        
        if fa == 0:
            s = a
            E = 0
            mensaje = str(s) + " es raiz de f(x)"
        
        elif fb == 0:
            s = b
            E = 0
            mensaje = str(s) + " es raiz de f(x)"
            
        elif fa*fb < 0:
            c = 0
            Xm = (a + b)/2                
            fe = func(Xm)
            fm.append(fe)
            x_vals.append(Xm)
            E.append(100)
            iters.append(c)
        
            while E[c] >= Tol and fe!= 0 and c < Niter:
                if fa*fe < 0: 
                    b = Xm              
                    fb = func(b)
                else:
                    a = Xm
                    fa = func(a)
                
                c = c+1
                Xm = (a + b)/2
                fe = func(Xm)
                fm.append(fe)
                x_vals.append(Xm)
                iters.append(c)
                
                if tipo_error == 1: 
                    Error_abs = abs(x_vals[c] - x_vals[c-1])
                    E.append(Error_abs)
                elif tipo_error == 2: 
                    Error_rel = abs((x_vals[c] - x_vals[c-1])/x_vals[c])
                    E.append(Error_rel)
                
            if fe == 0:
                s = Xm
                mensaje = str(s) + " es raiz de f(x)"
            
            elif E[c] < Tol:
                s = Xm
                d = {"Iteraciones": iters, "Xn": x_vals, "f(Xn)": fm, "Error": E}
                df = pd.DataFrame(d)

                mensaje = f"La solución aproximada es: {s}, con una tolerancia = {Tol} ({tipo})"
                
            else:
                s = Xm
                mensaje = "Fracasó en "+ str(Niter)+ " iteraciones "
        
        else:
            mensaje = "El intervalo es inadecuado"
            s=0


        graph = plot(Fun, (x,s-20,s+20) , xlabel='x', ylabel='y', show=False)
        graph.save('Grafica.png')

        context={'mensaje':mensaje,'df':df,'fun':Fun,'a':inia,'b':inib,'tipo_error':tipo,'num_tol':num_tol,'niter':Niter,'grafica':'../Grafica.png'}
        return render(request,template_name='1-biseccion.html',context=context)
    
    if request.method == 'GET':
        return render(request,template_name='1-biseccion.html')

def regla_falsa_page(request):
    x, y, z = symbols('x y z')
    if request.method == 'POST':
        Fun = request.POST.get('funcion')
        a = float(request.POST.get('a'))
        b = float(request.POST.get('b'))
        tipo_error = int(request.POST.get('tipo_de_tolerancia'))
        num_tol = float(request.POST.get('numero_de_tolerancia'))
        Niter = int(request.POST.get('iteraciones'))
        inia = a
        inib = b
        df=0

        if tipo_error == 1: 
            tipo = "decimales correctos"
            Tol = 0.5*(10**-num_tol)
        elif tipo_error == 2: 
            tipo = "cifras significativas"
            Tol = 5*(10**-num_tol)
        
        if a > b:
            temp = a
            a, b = b, temp
        
        Fun = sympify(Fun)
        func = lambda m: Fun.evalf(15, subs = {x: m})
        dom = continuous_domain(Fun, x, S.Reals)
        interval = Interval(a,b)
        if not interval.is_subset(dom):
            mensaje = f"La función no es continua en el intervalo [{a},{b}]. El método falla."
            graph = plot(Fun, xlabel='x', ylabel='y', show=False)
            graph.save('Grafica.png')
            context={'mensaje':mensaje,'df':df,'fun':Fun,'a':inia,'b':inib,'tipo_error':tipo,'num_tol':num_tol,'niter':Niter,'grafica':'../Grafica.png'}
            return render(request,template_name='2-regla_falsa.html',context=context)
        
        x_vals, fm, E, iters = [], [], [], []
        fa = func(a)
        fb = func(b)
        
        if fa == 0:
            s = a
            E = 0
            mensaje = str(s) + " es raiz de f(x)"
        
        elif fb == 0:
            s = b
            E = 0
            mensaje = str(s) + " es raiz de f(x)"
            
        elif fa*fb < 0:
            c = 0
            Xm = a - (fa*(b-a))/(fb - fa)         
            fe = func(Xm)
            fm.append(fe)
            x_vals.append(Xm)
            E.append(100)
            iters.append(c)
        
            while E[c] >= Tol and fe!= 0 and c < Niter:
                if fa*fe < 0: 
                    b = Xm              
                    fb = func(b)
                else:
                    a = Xm
                    fa = func(a)
                
                c = c+1
                Xm = a - (fa*(b-a))/(fb - fa)   
                fe = func(Xm)
                fm.append(fe)
                x_vals.append(Xm)
                iters.append(c)
                
                if tipo_error == 1: 
                    Error_abs = abs(x_vals[c] - x_vals[c-1])
                    E.append(Error_abs)
                elif tipo_error == 2: 
                    Error_rel = abs((x_vals[c] - x_vals[c-1])/x_vals[c])
                    E.append(Error_rel)
                
            if fe == 0:
                s = Xm
                mensaje = str(s) + " es raiz de f(x)"
            
            elif E[c] < Tol:
                s = Xm
                d = {"Iteraciones": iters, "Xn": x_vals, "f(Xn)": fm, "Error": E}
                df = pd.DataFrame(d)
                mensaje = f"La solución aproximada es: {s}, con una tolerancia = {Tol} ({tipo})"
           
            else:
                s = Xm
                mensaje = "Fracasó en "+ str(Niter)+ " iteraciones "
        
        else:
            mensaje = "El intervalo es inadecuado"
        
        graph = plot(Fun, (x,s-20,s+20), xlabel='x', ylabel='y', show=False)
        graph.save('Grafica.png')
        
        context={'mensaje':mensaje,'grafica':'../Grafica.png','df':df,'fun':Fun,'a':inia,'b':inib,'tipo_error':tipo,'num_tol':num_tol,'niter':Niter}
        return render(request,template_name='2-regla_falsa.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='2-regla_falsa.html')

def punto_fijo_page(request):
    x, y, z = symbols('x y z')
    if request.method == 'POST':
        Fun = request.POST.get('funcion')
        gf = request.POST.get('funcion_g')
        X0 = float(request.POST.get('x0'))
        tipo_error = int(request.POST.get('tipo_de_tolerancia'))
        num_tol = int(request.POST.get('numero_de_tolerancia'))
        Niter = int(request.POST.get('iteraciones'))
        iniX0 = X0
        df=0

        if tipo_error == 1: 
            tipo = "decimales correctos"
            Tol = 0.5*(10**-num_tol)
        elif tipo_error == 2: 
            tipo = "cifras significativas"
            Tol = 5*(10**-num_tol)
        
        Fun = sympify(Fun)
        gf = sympify(gf)
        func = lambda m: Fun.evalf(15, subs = {x: m})
        g = lambda m: gf.evalf(15, subs = {x: m})

        fn, xn, E, iters = [], [], [], []
        x_val = X0
        fe = func(X0)
        c, Error = 0, 100
        fn.append(fe)
        xn.append(X0)
        E.append(Error)
        iters.append(c)
        
        while E[c] >= Tol and fe!=0 and c < Niter:
            c = c+1
            x_val = g(X0)
            fe = func(x_val)
            fn.append(fe)
            xn.append(x_val)
            iters.append(c)
            
            if tipo_error == 1: 
                Error_abs = abs(xn[c] - xn[c-1])
                E.append(Error_abs)
            elif tipo_error == 2: 
                Error_rel = abs((xn[c] - xn[c-1])/xn[c])
                E.append(Error_rel)
        
            X0 = x_val
        
        if fe == 0:
            sol = x_val
            mensaje = str(sol) + " es raiz de f(x)"
            
        elif E[c] < Tol:
            sol = x_val 
            d = {"Iteraciones": iters, "Xn": xn, "f(Xn)": fn, "Error": E}
            df = pd.DataFrame(d)
            mensaje = f"La solución aproximada es: {sol}, con una tolerancia = {Tol} ({tipo})"
        
        else:
            sol = x_val
            mensaje = "Fracasó en "+ str(Niter)+ " iteraciones "
        
        graph = plot(Fun, (x,sol-20,sol+20), xlabel='x', ylabel='y', show=False)
        graph.save('Grafica.png')
        context={'mensaje':mensaje,'grafica':'../Grafica.png','df':df,'fun':Fun,'x0':iniX0,'gf':gf,'tipo_error':tipo,'num_tol':num_tol,'niter':Niter}
        return render(request,template_name='3-punto_fijo.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='3-punto_fijo.html')

def newton_page(request):
    x, y, z = symbols('x y z')
    if request.method == 'POST':
        Fun = request.POST.get('funcion')
        X0 = float(request.POST.get('x0'))
        tipo_error = int(request.POST.get('tipo_de_tolerancia'))
        num_tol = float(request.POST.get('numero_de_tolerancia'))
        Niter = int(request.POST.get('iteraciones'))
        iniX0 = X0
        df=0
        
        if tipo_error == 1: 
            tipo = "decimales correctos"
            Tol = 0.5*(10**-num_tol)
        elif tipo_error == 2: 
            tipo = "cifras significativas"
            Tol = 5*(10**-num_tol)

        Fun = sympify(Fun)
        func = lambda a: Fun.evalf(15, subs = {x: a})
        derivative = diff(Fun, x)
        deriv = lambda a: diff(Fun, x).evalf(15, subs = {x: a})

        domf = continuous_domain(Fun, x, S.Reals)
        print(domf)
        if not domf.contains(X0):
            mensaje = f"La función no está definida en x = {X0}. El método falla."
            graph = plot(Fun, xlabel='x', ylabel='y', show=False)
            graph.save('Grafica.png')
            context = {'mensaje':mensaje,'grafica':'../Grafica.png','df':df}
            return render(request,template_name='4-newton.html',context=context)
        
        elif not derivative.subs(x, X0).is_finite:
            mensaje = f"La función no es diferenciable en x = {X0}. El método falla."
            graph = plot(Fun, xlabel='x', ylabel='y', show=False)
            graph.save('Grafica.png')
            context = {'mensaje':mensaje,'grafica':'../Grafica.png','df':df}
            return render(request,template_name='4-newton.html',context=context)
        
        else: 
            fn, x_vals, E, dvs, iters = [], [], [], [], []
            xn = X0
            f, derivada = func(xn), deriv(xn)
            c, Error = 0, 100         
            fn.append(f)
            dvs.append(derivada)
            x_vals.append(xn)
            E.append(Error)
            iters.append(c)
            while (E[c] >= Tol) and (f != 0) and (derivada != 0) and (c < Niter):
                c = c+1
                xn = xn - (f/derivada)
                
                if not domf.contains(xn):
                    mensaje = f"La función no está definida en x{c} = {xn}. El método falla."
                    graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                    graph.save('Grafica.png')
                    context = {'mensaje':mensaje,'grafica':'../Grafica.png','df':df}
                    return render(request,template_name='4-newton.html',context=context)

                elif not derivative.subs(x,xn).is_finite:
                    mensaje = f"La función no es diferenciable en x{c} = {xn}. El método falla."
                    graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                    graph.save('Grafica.png')
                    context = {'mensaje':mensaje,'grafica':'../Grafica.png','df':df}
                    return render(request,template_name='4-newton.html',context=context)
                    
                f, derivada = func(xn), deriv(xn)
                fn.append(f)
                dvs.append(derivada)
                x_vals.append(xn)
                iters.append(c)
            
                if tipo_error == 1: 
                    Error_abs = abs(x_vals[c] - x_vals[c-1])
                    E.append(Error_abs)
                elif tipo_error == 2: 
                    Error_rel = abs((x_vals[c] - x_vals[c-1])/x_vals[c])
                    E.append(Error_rel)
            
            if f == 0:
                s = xn
                mensaje = str(s)+"es raiz de f(x)"
                graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                graph.save('Grafica.png')
                context = {'mensaje':mensaje,'grafica':'../Grafica.png','df':df}
                return render(request,template_name='4-newton.html',context=context)
            
            if derivada == 0:
                s = xn
                mensaje = f"La derivada es 0 en x = {s}. El método falla."
                graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                graph.save('Grafica.png')
                context = {'mensaje':mensaje,'grafica':'../Grafica.png','df':df}
                return render(request,template_name='4-newton.html',context=context)
                
            elif E[c] < Tol:
                s = xn
                d = {"Iteraciones": iters, "Xn": x_vals, "f(Xn)": fn, "f'(Xn)": dvs, "Error": E}
                df = pd.DataFrame(d) 
                mensaje = f"La solución aproximada es: {s}, con una tolerancia = {Tol} ({tipo})"
                graph = plot(Fun, (x,s-20,s+20), xlabel='x', ylabel='y', show=False)
                graph.save('Grafica.png')
                context={'mensaje':mensaje,'grafica':'../Grafica.png','df':df,'fun':Fun,'x0':iniX0,'tipo_error':tipo,'num_tol':num_tol,'niter':Niter}
                return render(request,template_name='4-newton.html',context=context)
            
            else:
                s = xn
                mensaje = "Fracasó en"+ str(Niter) + "iteraciones"
                graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                graph.save('Grafica.png')
                context = {'mensaje':mensaje,'grafica':'../Grafica.png','df':df}
                return render(request,template_name='4-newton.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='4-newton.html')
    
def raices_multiples_page(request):
    x, y, z = symbols('x y z')
    if request.method == 'POST':
        Fun = request.POST.get('funcion')
        X0 = float(request.POST.get('x0'))
        tipo_error = int(request.POST.get('tipo_de_tolerancia'))
        num_tol = float(request.POST.get('numero_de_tolerancia'))
        Niter = int(request.POST.get('iteraciones'))
        iniX0 = X0
        df=0
        
        if tipo_error == 1: 
            tipo = "decimales correctos"
            Tol = 0.5*(10**-num_tol)
        elif tipo_error == 2: 
            tipo = "cifras significativas"
            Tol = 5*(10**-num_tol)
        
        Fun = sympify(Fun)
        func = lambda a: Fun.evalf(15, subs = {x: a})
        derivative = diff(Fun, x)
        deriv = lambda a: diff(Fun, x).evalf(15, subs = {x: a})
        derivative2 = diff(Fun, x, 2)
        deriv2 = lambda a: diff(Fun, x, 2).evalf(15, subs = {x: a})
        
        domf = continuous_domain(Fun, x, S.Reals)
        if not domf.contains(X0):
            mensaje = f"La función no está definida en x = {X0}. El método falla."
            graph = plot(Fun, xlabel='x', ylabel='y', show=False)
            graph.save('Grafica.png')
            context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
            return render(request,template_name='5-raices_multiples.html',context=context)
            
        elif not derivative.subs(x, X0).is_finite:
            mensaje = f"La función no es diferenciable en x = {X0}. El método falla."
            graph = plot(Fun, xlabel='x', ylabel='y', show=False)
            graph.save('Grafica.png')
            context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
            return render(request,template_name='5-raices_multiples.html',context=context)
        
        elif not derivative2.subs(x, X0).is_finite:
            mensaje = f"La función no tiene segunda derivada en x = {X0}. El método falla."
            graph = plot(Fun, xlabel='x', ylabel='y', show=False)
            graph.save('Grafica.png')
            context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
            return render(request,template_name='5-raices_multiples.html',context=context)
        
        else: 
            fn, x_vals, E, dvs, dvs2, iters = [], [], [], [], [], []
            xn = X0
            f, derivada, derivada2 = func(xn), deriv(xn), deriv2(xn)
            c, Error = 0, 100         
            fn.append(f)
            dvs.append(derivada)
            dvs2.append(derivada2)
            x_vals.append(xn)
            E.append(Error)
            iters.append(c)
            
            while (E[c] >= Tol) and (f != 0) and (derivada**2 - (f*derivada2)) != 0  and (c < Niter):
                c = c + 1
                xn = xn - (f*derivada)/(derivada**2 - (f*derivada2))
                
                if not domf.contains(xn):
                    mensaje = f"La función no está definida en x{c} = {xn}. El método falla."
                    graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                    graph.save('Grafica.png')
                    context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                    return render(request,template_name='5-raices_multiples.html',context=context)

                elif not derivative.subs(x,xn).is_finite:
                    mensaje = f"La función no es diferenciable en x{c} = {xn}. El método falla."
                    graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                    graph.save('Grafica.png')
                    context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                    return render(request,template_name='5-raices_multiples.html',context=context)
                
                elif not derivative2.subs(x,xn).is_finite:
                    mensaje = f"La función no tiene segunda derivada en x{c} = {xn}. El método falla."
                    graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                    graph.save('Grafica.png')
                    context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                    return render(request,template_name='5-raices_multiples.html',context=context)
                
                f, derivada, derivada2 = func(xn), deriv(xn), deriv2(xn)
                fn.append(f)
                dvs.append(derivada)
                dvs2.append(derivada2)
                x_vals.append(xn)
                iters.append(c)
                
                if tipo_error == 1: 
                    Error_abs = abs(x_vals[c] - x_vals[c-1])
                    E.append(Error_abs)
                elif tipo_error == 2: 
                    Error_rel = abs((x_vals[c] - x_vals[c-1])/x_vals[c])
                    E.append(Error_rel)
            
            if f == 0:
                s = xn
                mensaje = str(s)+"es raiz de f(x)"
                graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                graph.save('Grafica.png')
                context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                return render(request,template_name='5-raices_multiples.html',context=context)
            
            if (derivada**2 - (f*derivada2)) == 0:
                s = xn
                mensaje = "El método falla. El denominador es 0."
                graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                graph.save('Grafica.png')
                context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                return render(request,template_name='5-raices_multiples.html',context=context)
                
            elif E[c] < Tol: #Vuelve
                s = xn
                d = {"Iteraciones": iters, "Xn": x_vals, "f(Xn)": fn, "f'(Xn)": dvs, "f''(Xn)": dvs2, "Error": E}
                df = pd.DataFrame(d) 
                mensaje = f"La solución aproximada es: {s}, con una tolerancia = {Tol} ({tipo})"
                graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                graph.save('Grafica.png')
                context={'mensaje':mensaje,'grafica':'../Grafica.png','df':df,'fun':Fun,'x0':iniX0,'tipo_error':tipo,'num_tol':num_tol,'niter':Niter}
                return render(request,template_name='5-raices_multiples.html',context=context)
            
            else:
                s = xn
                mensaje = f"Fracasó en {Niter} iteraciones"
                graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                graph.save('Grafica.png')
                context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                return render(request,template_name='5-raices_multiples.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='5-raices_multiples.html')

def secante_page(request):
    if request.method == 'POST':
        x, y, z = symbols('x y z')
        Fun = request.POST.get('funcion')
        X0 = float(request.POST.get('x0'))
        X1 = float(request.POST.get('x1'))
        tipo_error = int(request.POST.get('tipo_de_tolerancia'))
        num_tol = float(request.POST.get('numero_de_tolerancia'))
        Niter = int(request.POST.get('iteraciones'))
        iniX0 = X0
        iniX1 = X1
        df=0
        
        if tipo_error == 1: 
            tipo = "decimales correctos"
            Tol = 0.5*(10**-num_tol)
        elif tipo_error == 2: 
            tipo = "cifras significativas"
            Tol = 5*(10**-num_tol)
        
        Fun = sympify(Fun)
        func = lambda m: Fun.evalf(15, subs = {x: m})
        
        domf = continuous_domain(Fun, x, S.Reals)
        if not domf.contains(X0):
            mensaje = f"La función no está definida en x = {X0}. El método falla."
            graph = plot(Fun, xlabel='x', ylabel='y', show=False)
            graph.save('Grafica.png')
            context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
            return render(request,template_name='6-secante.html',context=context)
        elif not domf.contains(X1):
            mensaje = f"La función no está definida en x = {X1}. El método falla."
            graph = plot(Fun, xlabel='x', ylabel='y', show=False)
            graph.save('Grafica.png')
            context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
            return render(request,template_name='6-secante.html',context=context)
        
        else: 
            fn, xn, E, iters = [], [], [], []
            f0, f1 = func(X0), func(X1)
            fn.extend([f0,f1])
            xn.extend([X0,X1])
            E.extend([abs(100), abs(X0-X1)])
            iters.extend([0,1])
            
            if f0 == 0:
                mensaje = str(X0)+" es raiz de f(x)"
            
            elif f1 == 0:
                mensaje = str(X1)+" es raiz de f(x)"
            
            else:
                c = 1
                x_val = X1
                while (E[c] >= Tol) and (fn[c] != 0) and (f1-f0 != 0)  and (c < Niter):
                    X1, X0 = xn[c], xn[c-1]
                    f1, f0 = fn[c], fn[c-1]
                    c = c+1
                    x_val = X1 - (f1*(X1-X0))/(f1 - f0)   
                    
                    if not domf.contains(x_val):
                        mensaje = f"La función no está definida en x{c} = {x_val}. El método falla."
                        graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                        graph.save('Grafica.png')
                        context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                        return render(request,template_name='6-secante.html',context=context)
                        
                        
                    fe = func(x_val)
                    fn.append(fe)
                    xn.append(x_val)
                    iters.append(c)
        
                    if tipo_error == 1: 
                        Error_abs = abs(xn[c] - xn[c-1])
                        E.append(Error_abs)
                    elif tipo_error == 2: 
                        Error_rel = abs((xn[c] - xn[c-1])/xn[c])
                        E.append(Error_rel)
                
                if fn[c] == 0:
                    s = x_val
                    mensaje = str(s)+"es raiz de f(x)"
                    graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                    graph.save('Grafica.png')
                    context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                    return render(request,template_name='6-secante.html',context=context)
                
                if (f1-f0 == 0):
                    mensaje = "El método falla"
                    graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                    graph.save('Grafica.png')
                    context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                    return render(request,template_name='6-secante.html',context=context)
                
                elif E[c] < Tol:
                    s = x_val
                    d = {"Iteraciones": iters, "Xn": xn, "f(Xn)": fn, "Error": E}
                    df = pd.DataFrame(d) 
                    mensaje = f"La solución aproximada es: {s}, con una tolerancia = {Tol} ({tipo})"
                    graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                    graph.save('Grafica.png')
                    context={'mensaje':mensaje,'grafica':'../Grafica.png','df':df,'fun':Fun,'x0':iniX0,'x1':iniX1,'tipo_error':tipo,'num_tol':num_tol,'niter':Niter}
                    return render(request,template_name='6-secante.html',context=context)
                
                else:
                    s = x_val
                    mensaje = f"Fracasó en {Niter} iteraciones "
                    graph = plot(Fun, xlabel='x', ylabel='y', show=False)
                    graph.save('Grafica.png')
                    context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png'}
                    return render(request,template_name='6-secante.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='6-secante.html')

def jacobi_page(request):

    if request.method == 'POST':
        size = int(request.POST.get('tamaño'))
        input_data = request.POST.get('matriz')
        b = request.POST.get('b')
        x0 = request.POST.get('x0')
        tipo_error = int(request.POST.get('tipo_de_tolerancia')) #1 = decimales correctos, 2 = cifras significativas
        num_tol = float(request.POST.get('numero_de_tolerancia')) #numero de decimales correctos o cifras significativas
        niter = int(request.POST.get('iteraciones'))
        norma = request.POST.get('norma')  #Escribir que puede ser un numero o inf 
        df=0
        iniX0=x0
        
        #---------- Creacion de A
        b = np.asarray(list(map(float, b.split()))).T
        x0 = np.asarray(list(map(float, x0.split()))).T
        
        if len(b) != size : 
            mensaje = f"Error: El vector b debe tener un tamaño de {size} componentes"
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='7-jacobi.html',context=context)
        if len(x0) != size : 
            mensaje = f"Error: El vector inicial (x0) debe tener un tamaño de {size} componentes"
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='7-jacobi.html',context=context)
        
        rows = input_data.split(';')
        A, A_flag=[], True
        for row in rows:
            row_data = row.strip().split()
            row_numbers = list(map(float, row_data))
            if len(row_numbers) != size: A_flag= False
            A.append(row_numbers)
        if len(A) != size: A_flag = False
        A = np.array(A)
        if not A_flag: 
            mensaje = f"Error: La matriz A no corresponde una una matriz cuadrada de dimensión {size} "
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='7-jacobi.html',context=context)
    # -----------------------------------------

        if tipo_error == 1: 
            tipo = "decimales correctos"
            Tol = 0.5*(10**-num_tol)
        elif tipo_error == 2: 
            tipo = "cifras significativas"
            Tol = 5*(10**-num_tol)
        
            
        row_sums = np.sum(np.abs(A), axis=1)
        diag = np.abs(np.diag(A))
        dominant = np.all(diag > row_sums - diag)
        
        if norma =="inf": norma=np.inf
        elif int(norma) >=1: norma = int(norma)
        else : norma = np.inf
        
        c = 0
        error = Tol + 1
        D = np.diag(np.diag(A))
        L = -np.tril(A, -1)
        U = -np.triu(A, 1)
        E, xn = [],[]
        
        if(np.linalg.det(D)==0) : 
            mensaje = "Matriz D es singular, el método falla" ### error
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='7-jacobi.html',context=context)
        else:
            T = np.linalg.inv(D) @ (L + U)
            C = np.linalg.inv(D) @ b
            while error > Tol and c < niter:
                x1 = T @ x0 + C
                error = np.linalg.norm(x1 - x0, norma)
                if tipo_error == 2:
                    error = error/np.linalg.norm(x1, norma)
                E.append(error)
                xn.append(x0)
                x0 = x1
                c += 1
            
            eigenvals = np.linalg.eigvals(T)
            rho = np.max(np.abs(eigenvals))   
            if error < Tol:
                s = x0
                error = np.linalg.norm(x1 - x0, norma)
                if tipo_error == 2:
                    error = error/np.linalg.norm(x1, norma)
                E.append(error)
                xn.append(x0)
                itera = np.arange(0,c+1,1)
                df = pd.DataFrame({"Iteraciones": itera,"Xn": xn, "Error": E})     
                mensaje = f"La solución aproximada es: {np.ravel(s)}, con una tolerancia = {Tol} ({tipo})"
                if rho < 1: mensaje = mensaje + "\n "+ f"Esta solucion es única porque el radio espectral de T es {rho} y es menor a 1"
                context={'mensaje':mensaje,'df':df,'tamaño':size,'matriz':input_data,'x0':iniX0,'b':b,'tipo_error':tipo,'num_tol':num_tol,'niter':niter,'norma':norma}
                return render(request,template_name='7-jacobi.html',context=context)
            else:
                s = x0
                print(f"Fracasó en {niter} iteraciones")
                if rho >= 1:
                    mensaje = f"Es posible que haya fallado el método porque el radio espectral de T es {rho} y es mayor a 1"
                    context = {'mensaje':mensaje,'df':df}
                    return render(request,template_name='7-jacobi.html',context=context)
                    
                if not dominant:
                    mensaje = "Es posible que haya fallado el método porque A no es diagonal dominante"
                    context = {'mensaje':mensaje,'df':df}
                    return render(request,template_name='7-jacobi.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='7-jacobi.html')

def gauss_seidel_page(request):
    if request.method == 'POST':
        size = int(request.POST.get('tamaño'))
        input_data = request.POST.get('matriz')
        b = request.POST.get('b')
        x0 = request.POST.get('x0')
        tipo_error = int(request.POST.get('tipo_de_tolerancia')) #1 = decimales correctos, 2 = cifras significativas
        num_tol = float(request.POST.get('numero_de_tolerancia')) #numero de decimales correctos o cifras significativas
        niter = int(request.POST.get('iteraciones'))
        norma = request.POST.get('norma')  #Escribir que puede ser un numero o inf 
        df=0
        iniX0=x0

        b = np.asarray(list(map(float, b.split()))).T
        x0 = np.asarray(list(map(float, x0.split()))).T
        
        if len(b) != size :
            mensaje = f"Error: El vector b debe tener un tamaño de {size} componentes"
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='8-gauss_seidel.html',context=context)
        if len(x0) != size : 
            mensaje = f"Error: El vector inicial (x0) debe tener un tamaño de {size} componentes"
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='8-gauss_seidel.html',context=context)
        rows = input_data.split(';')
        A, A_flag=[], True
        for row in rows:
            row_data = row.strip().split()
            row_numbers = list(map(float, row_data))
            if len(row_numbers) != size: A_flag= False
            A.append(row_numbers)

        if len(A) != size: A_flag = False
        A = np.array(A)
        if not A_flag: 
            mensaje = f"Error: La matriz A no corresponde a una matriz cuadrada de dimensión {size} "
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='8-gauss_seidel.html',context=context)
                
        if tipo_error == 1: 
            tipo = "decimales correctos"
            Tol = 0.5*(10**-num_tol)
        elif tipo_error == 2: 
            tipo = "cifras significativas"
            Tol = 5*(10**-num_tol)
        
        if norma =="inf": norma=np.inf
        elif int(norma) >=1: norma = int(norma)
        else : norma = np.inf
        
        row_sums = np.sum(np.abs(A), axis=1)
        diag = np.abs(np.diag(A))
        dominant = np.all(diag > row_sums - diag)

        c = 0
        error = Tol + 1
        D = np.diag(np.diag(A))
        L = -np.tril(A, -1)
        U = -np.triu(A, 1)
        E = np.zeros(niter)
        
        
        x_vals, E, iters = [], [], []
        x_vals.append(np.ravel(x0))
        E.append(error)
        iters.append(c)
        
        if np.linalg.det(D - L) == 0: 
            mensaje = "La matriz (D - L) no es invertible. El método falla"
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='8-gauss_seidel.html',context=context)
        else:
            T = np.linalg.inv(D - L) @ U
            C = np.linalg.inv(D - L) @ b
            eigenvalues = np.linalg.eigvals(T)
            radio_espectral = max(abs(eigenvalues))
            while error > Tol and c < niter:
                x1 = T @ x0 + C
                x_vals.append(np.ravel(x1))
                error = np.linalg.norm(x1 - x0, norma)
                
                if tipo_error == 2:
                    error = error/np.linalg.norm(x1, norma)
            
                E.append(error)
                x0 = x1
                c += 1
                iters.append(c)

            if error < Tol:
                s = x0
                d = {"Iteraciones": iters, "Xn": x_vals, "Error": E}
                df = pd.DataFrame(d)  
                mensaje = f"La solución aproximada es: {np.ravel(s)}, con una tolerancia = {Tol} ({tipo})"
                if radio_espectral < 1:
                    mensaje = mensaje + "\n "+ f"Esta solucion es única porque el radio espectral de T es {radio_espectral} y es menor a 1"
                context={'mensaje':mensaje,'df':df,'tamaño':size,'matriz':input_data,'x0':iniX0,'b':b,'tipo_error':tipo,'num_tol':num_tol,'niter':niter,'norma':norma}
                return render(request,template_name='8-gauss_seidel.html',context=context)
            else:
                s = x0
                print(f'Fracasó en {niter} iteraciones')
                if radio_espectral >= 1:
                    mensaje = f"Es posible que haya fallado el método porque el radio espectral de T es {radio_espectral} y es mayor a 1"
                    context = {'mensaje':mensaje,'df':df}
                    return render(request,template_name='8-gauss_seidel.html',context=context)

                if not dominant: 
                    mensaje = "Es posible que haya fallado el método porque A no es diagonal dominante"
                    context = {'mensaje':mensaje,'df':df}
                    return render(request,template_name='8-gauss_seidel.html',context=context)


    if request.method == 'GET':
        return render(request,template_name='8-gauss_seidel.html')

def SOR_page(request):
    if request.method == 'POST':
        size = int(request.POST.get('tamaño'))
        input_data = request.POST.get('matriz')
        b = request.POST.get('b')
        x0 = request.POST.get('x0')
        w = float(request.POST.get('w'))
        tipo_error = int(request.POST.get('tipo_de_tolerancia')) #1 = decimales correctos, 2 = cifras significativas
        num_tol = float(request.POST.get('numero_de_tolerancia')) #numero de decimales correctos o cifras significativas
        niter = int(request.POST.get('iteraciones'))
        norma = request.POST.get('norma')  #Escribir que puede ser un numero o inf 
        df=0
        iniX0=x0

        b = np.asarray(list(map(float, b.split()))).T
        x0 = np.asarray(list(map(float, x0.split()))).T
        
        if len(b) != size : 
            mensaje = f"Error: El vector b debe tener un tamaño de {size} componentes"
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='9-SOR.html',context=context)
        if len(x0) != size : 
            mensaje = f"Error: El vector inicial (x0) debe tener un tamaño de {size} componentes"
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='9-SOR.html',context=context)
        rows = input_data.split(';')
        A, A_flag=[], True
        for row in rows:
            row_data = row.strip().split()
            row_numbers = list(map(float, row_data))
            if len(row_numbers) != size: A_flag= False
            A.append(row_numbers)
        
        if len(A) != size: A_flag = False
        A = np.array(A)
        if not A_flag: 
            mensaje = f"Error: La matriz A no corresponde a una matriz cuadrada de dimensión {size} "
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='9-SOR.html',context=context)
        
        if w<=0 or w>=2: 
            mensaje = " El peso de w no está dentro de su dominio. Elige un número entre I = (0,2) para que el método converja" #Imprimir este error
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='9-SOR.html',context=context)

        if tipo_error == 1: 
            tipo = "decimales correctos"
            Tol = 0.5*(10**-num_tol)
        elif tipo_error == 2: 
            tipo = "cifras significativas"
            Tol = 5*(10**-num_tol)
            
        if norma =="inf": norma=np.inf
        elif int(norma) >=1: norma = int(norma)
        else : norma = np.inf
        
        row_sums = np.sum(np.abs(A), axis=1)
        diag = np.abs(np.diag(A))
        dominant = np.all(diag > row_sums - diag)
        
        simetria = np.all(np.abs(A-A.T) < 1e-8)
        def_positiva = np.all(np.linalg.eigvals(A) > 0)
        
        c = 0
        error = Tol + 1
        D = np.diag(np.diag(A))
        L = -np.tril(A, k=-1)
        U = -np.triu(A, k=1)
        E = [Tol + 1]
        x_n = [x0]
        
        if np.linalg.det(D - w * L)== 0:  
            mensaje = "Matriz D-wL es singular, el método falla"
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='9-SOR.html',context=context)
        
        T = np.linalg.inv(D - w * L) @ ((1 - w) * D + w * U)
        C = w * np.linalg.inv(D - w * L) @ b
        
        eigenvals = np.linalg.eigvals(T)
        rho = np.max(np.abs(eigenvals))
        while error > Tol and c < niter:
            x1 = T @ x0 + C
            error = np.linalg.norm(x1 - x0, norma)
            if tipo_error == 2:
                error = error/np.linalg.norm(x1, norma)
            E.append(error)
            x0 = x1
            x_n.append(x0)
            c += 1
        if error < Tol:
            s = x0
            itera = np.arange(0,c+1,1)
            df = pd.DataFrame({"Iteraciones": itera,"Xn": x_n, "Error":E})
            print(df)
            mensaje = f'La solución aproximada es: {np.ravel(s)}, con una tolerancia = {Tol} ({tipo})'
            if rho < 1: mensaje = mensaje + f"Esta solucion es única porque el radio espectral de T es {rho} y es menor a 1"
            context={'mensaje':mensaje,'df':df,'tamaño':size,'matriz':input_data,'x0':iniX0,'b':b,'tipo_error':tipo,'num_tol':num_tol,'niter':niter,'w':w,'norma':norma}
            return render(request,template_name='9-SOR.html',context=context)
        else:
            s = x0
            mensaje = f'Fracasó en {niter} iteraciones. '
            if not dominant:
                mensaje = mensaje + "El método puede diverger porque la matriz A no es estrictamente diagonal dominante. \n"
            mensaje = mensaje + "Es posible que haya fallado el método porque: \n"
            if rho >= 1:
                mensaje = mensaje + f"- El radio espectral de T es {rho} y es mayor a 1. "
            if not simetria: 
                mensaje = mensaje + "- La matriz A no es simétrica. "
            if not def_positiva:
                mensaje = mensaje + "- La matriz A no es definida positiva. "
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='9-SOR.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='9-SOR.html')

def vandermonde_page(request):
    if request.method == 'POST':
        x = request.POST.get('x')
        y = request.POST.get('y')
        
        x = np.asarray(list(map(float, x.split())))
        y = np.asarray(list(map(float, y.split())))

        if len(x) != len(y):
            mensaje = "Los valores de X y Y deben tener la misma cantidad de elementos "
            context = {'mensaje':mensaje}
            return render(request,template_name='10-vandermonde.html',context=context)

        
        if len(np.unique(x)) == len(x): different = True
        else: different = False
        if not different: 
            mensaje = "Los valores de x no deben repertirse para que sea una función "
            context = {'mensaje':mensaje}
            return render(request,template_name='10-vandermonde.html',context=context)

        
        n = len(x)
        A = np.zeros((n, n))
        sorted_indices = np.argsort(x)
        x = x[sorted_indices]
        y = y[sorted_indices]

        for i in range(n):
            A[:,i] = np.power(x, n-i-1)
        b = y.T
        if np.linalg.det(A)== 0: 
            mensaje = "El método no funciona porque la matriz A es singular"
            context = {'mensaje':mensaje}
            return render(request,template_name='10-vandermonde.html',context=context)
        else:
            pol = np.linalg.solve(A, b)
            #Imprimir el polinomio---------------------------------------
            pol_string, i = str(round(pol[0],3)) + " x^" + str(len(pol)-1), 1
            for comp in pol[1:-1]:
                if round(comp,3)>=0: pol_string+= " + " + str(round(comp,3)) + " x^" + str(len(pol)-1-i) +" "
                else : pol_string+= str(round(comp,3)) + " x^" + str(len(pol)-1-i) +" "
                i+=1
            if round(pol[-1],3)>=0: pol_string+=" + "+str(round(pol[-1],3))
            else: pol_string+= str(round(pol[-1],3))
            #----------------------------------------------------------------
            mensaje = f"El polinomio que interpola los puntos dados es: {pol_string}"
            #Grafica del polinomio ------------------------------------------
            z, paso = 0, 0.1
            t = np.arange(min(x),max(x)+paso,paso) #domain function
            for i, val in zip(range(len(pol)), pol):
                z += val*t**(len(pol)-i-1)

            plt.cla()
            plt.clf()
            plt.plot(t,z)
            plt.plot(x,y, marker='.', ls='none', ms=10, color="k")
            plt.title("Gráfica de la interpolación de los puntos")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.legend(["Interpolación","Puntos dados"])
            plt.grid()
            plt.savefig('Grafica.png')
            plt.cla()
            plt.clf()

            context = {'mensaje':mensaje,'grafica':'../Grafica.png','x':x,'y':y,'paso':True}
            return render(request,template_name='10-vandermonde.html',context=context)
    
    if request.method == 'GET':
        return render(request,template_name='10-vandermonde.html')

def newton_diferencias_divididas_page(request):
    if request.method == 'POST':
        x = request.POST.get('x')
        y = request.POST.get('y')
        df=0
        
        x = np.asarray(list(map(float, x.split())))
        y = np.asarray(list(map(float, y.split())))

        if len(x) != len(y):
            mensaje = "Los valores de X y Y deben tener la misma cantidad de elementos "
            context = {'mensaje':mensaje}
            return render(request,template_name='11-newton_diferencias_divididas.html',context=context)
        
    
        if len(np.unique(x)) == len(x): different = True
        else: different = False
        if not different: 
            mensaje = "Los valores de x no deben repertirse para que sea una función " 
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='11-newton_diferencias_divididas.html',context=context)
    
        n = len(x)
        sorted_indices = np.argsort(x)
        x = x[sorted_indices]
        y = y[sorted_indices]
        Tabla = np.zeros((n, n + 1))
        Tabla[:, 0] = x
        Tabla[:, 1] = y
        
        for j in range(2, n + 1):
            for i in range(j - 1, n):
                Tabla[i, j] = (Tabla[i, j - 1] - Tabla[i - 1, j - 1]) / (Tabla[i, 0] - Tabla[i - j + 1, 0])
        df = pd.DataFrame(Tabla)
        coef = np.diag(Tabla[:, 1:])
        #######Newtonor
        pol = np.array([1])
        acum = pol
        pol = coef[0] * acum
            
        for i in range(n - 1):
            pol = np.concatenate(([0], pol))
            acum = np.convolve(acum, [1, -x[i]])
            pol = pol + coef[i + 1] * acum
        
        # Imprimir el polinomio ---------------------------------------------
        pol_string, i = str(round(pol[0],3)) + " x^" + str(len(pol)-1), 1
        for comp in pol[1:-1]:
            if round(comp,3)>=0: pol_string+= " + " + str(round(comp,3)) + " x^" + str(len(pol)-1-i) +" "
            else : pol_string+= str(round(comp,3)) + " x^" + str(len(pol)-1-i) +" "
            i+=1
        if round(pol[-1],3)>=0: pol_string+=" + "+str(round(pol[-1],3))
        else: pol_string+= str(round(pol[-1],3))
        
        mensaje = f"El polinomio que interpola los puntos dados es: {pol_string}"
        # Grafica -------------------------------------------------------------
            
        z, paso = 0, 0.1
        t = np.arange(min(x),max(x)+paso,paso) #domain function
        for i, val in zip(range(len(pol)), pol):
            z += val*t**(len(pol)-i-1)

        plt.cla()
        plt.clf()
        plt.plot(x,y, marker='.', ls='none', ms=10, color="k")
        plt.plot(t,z)
        plt.title(f"Gráfica de la interpolación de los puntos")
        plt.legend(["Puntos dados", "Interpolación"])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid()
        plt.savefig('Grafica.png')
        plt.cla()
        plt.clf()

        context = {'mensaje':mensaje,'grafica':'../Grafica.png','x':x,'y':y,'paso':True,'df':df}
        return render(request,template_name='11-newton_diferencias_divididas.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='11-newton_diferencias_divididas.html')

def lagrange_page(request):
    if request.method == 'POST':
        x = request.POST.get('x')
        y = request.POST.get('y')
        
        x = np.asarray(list(map(float, x.split())))
        y = np.asarray(list(map(float, y.split())))

        if len(x) != len(y):
            mensaje = "Los valores de X y Y deben tener la misma cantidad de elementos "
            context = {'mensaje':mensaje}
            return render(request,template_name='12-lagrange.html',context=context)
        
        if len(np.unique(x)) == len(x): different = True
        else: different = False
        if not different: 
            mensaje = "Los valores de x no deben repertirse para que sea una función "
            context = {'mensaje':mensaje}
            return render(request,template_name='12-lagrange.html',context=context)
        
        n = len(x)
        Tabla = np.zeros((n, n))
        sorted_indices = np.argsort(x)
        x = x[sorted_indices]
        y = y[sorted_indices]
        for i in range(n):
            Li = 1
            den = 1
            for j in range(n):
                if j != i:
                    paux = [1, -x[j]]
                    Li = np.convolve(Li, paux)
                    den *= (x[i] - x[j])
            Tabla[i, :] = y[i] * Li / den
        pol = np.sum(Tabla, axis=0)
        
        #Imprimir el polinomio---------------------------------------
        pol_string, i = str(round(pol[0],3)) + " x^" + str(len(pol)-1), 1
        for comp in pol[1:-1]:
            if round(comp,3)>=0: pol_string+= " + " + str(round(comp,3)) + " x^" + str(len(pol)-1-i) +" "
            else : pol_string+= str(round(comp,3)) + " x^" + str(len(pol)-1-i) +" "
            i+=1
        if round(pol[-1],3)>=0: pol_string+=" + "+str(round(pol[-1],3))
        else: pol_string+= str(round(pol[-1],3))
        #----------------------------------------------------------------
        mensaje = f"El polinomio que interpola los puntos dados es: {pol_string}"
        #Grafica del polinomio ------------------------------------------
        z, paso = 0, 0.1
        t = np.arange(min(x),max(x)+paso,paso) #domain function
        for i, val in zip(range(len(pol)), pol):
            z += val*t**(len(pol)-i-1)
        
        plt.cla()
        plt.clf()
        plt.plot(x,y, marker='.', ls='none', ms=10, color="k")
        plt.plot(t,z)
        plt.title(f"Gráfica de la interpolación de los puntos")
        plt.legend(["Puntos dados", "Interpolación"])
        plt.grid()
        plt.xlabel("x")
        plt.ylabel("y")
        plt.savefig('Grafica.png')
        plt.cla()
        plt.clf()

        context = {'mensaje':mensaje,'grafica':'../Grafica.png','x':x,'y':y,'paso':True}
        return render(request,template_name='12-lagrange.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='12-lagrange.html')

def spline_page(request):
    if request.method == 'POST':
        x = request.POST.get('x')
        y = request.POST.get('y')
        d = int(request.POST.get('tipo_de_spline')) # 1: lineal, 3: cubico
        df=0
        mensaje=""
        tipo=""

        if d==1:
            tipo="Lineal"
        elif d==3:
            tipo="Cubica"
        
        x = np.asarray(list(map(float, x.split())))
        y = np.asarray(list(map(float, y.split())))

        if len(x) != len(y):
            mensaje = "Los valores de X y Y deben tener la misma cantidad de elementos "
            context = {'mensaje':mensaje}
            return render(request,template_name='13-spline.html',context=context)
    
        if len(np.unique(x)) == len(x): different = True
        else: different = False
        if not different: 
            mensaje = "Los valores de x no deben repertirse para que sea una función " #error
            context = {'mensaje':mensaje,'df':df}
            return render(request,template_name='13-spline.html',context=context)

        n = len(x)
        
        # Ordenar los vectores
        sorted_indices = np.argsort(x)
        x = x[sorted_indices]
        y = y[sorted_indices]
        
        A = np.zeros(((d+1)*(n-1),(d+1)*(n-1)))
        b = np.zeros((d+1)*(n-1))
        
        cua = np.power(x, 2)
        cub = np.power(x, 3)
        flag, bad = False, False
        # Lineal
        if d == 1:
            c = 0
            h = 0
            for i in range(n-1):
                A[i, c] = x[i]
                A[i, c+1] = 1
                b[i] = y[i]
                c += 2
                h += 1
            
            c = 0
            for i in range(1, n):
                A[h, c] = x[i]
                A[h, c+1] = 1
                b[h] = y[i]
                c += 2
                h += 1
                
        # Cubic
        elif d == 3:
            c, h = 0, 0
            for i in range(n-1):
                A[i, c] = cub[i]
                A[i, c+1] = cua[i]
                A[i, c+2] = x[i]
                A[i, c+3] = 1
                b[i] = y[i]
                c += 4
                h += 1
            
            c = 0
            for i in range(1, n):
                A[h, c] = cub[i]
                A[h, c+1] = cua[i]
                A[h, c+2] = x[i]
                A[h, c+3] = 1
                b[h] = y[i]
                c += 4
                h += 1
            
            c = 0
            for i in range(1, n-1):
                A[h, c] = 3*cua[i]
                A[h, c+1] = 2*x[i]
                A[h, c+2] = 1
                A[h, c+4] = -3*cua[i]
                A[h, c+5] = -2*x[i]
                A[h, c+6] = -1
                b[h] = 0
                c += 4
                h += 1
            
            c = 0
            for i in range(1, n-1):
                A[h, c] = 6*x[i]
                A[h, c+1] = 2
                A[h, c+4] = -6*x[i]
                A[h, c+5] = -2
                b[h] = 0
                c += 4
                h += 1
            
            A[h, 0] = 6*x[0]
            A[h, 1] = 2
            b[h] = 0
            h += 1
            A[h, c] = 6*x[-1]
            A[h, c+1] = 2
            b[h] = 0
        if not flag:
            val = np.linalg.inv(A).dot(b)
            reshaped_matrix = val.reshape(len(x)-1, d+1)
            Tabla = reshaped_matrix
            #---Imprimir polinomio --------
            Tabla_string, j= [0]*len(Tabla),0
            l=1
            for pol in Tabla:
                pol_string, i = str(round(pol[0],3)) + " x^" + str(len(pol)-1), 1
                for comp in pol[1:-1]:
                    if round(comp,3)>=0: pol_string+= " + " + str(round(comp,3)) + " x^" + str(len(pol)-1-i) +" "
                    else : pol_string+= str(round(comp,3)) + " x^" + str(len(pol)-1-i) +" "
                    i+=1
                    
                l+=1
                if round(pol[-1],3)>=0: pol_string+=" + "+str(round(pol[-1],3))
                else: pol_string+= str(round(pol[-1],3))
                
                Tabla_string[j]= pol_string
                j+=1
            mensaje = mensaje + f"Los polinomios de grado {d} que interpolan los puntos son :"
            df = pd.DataFrame({"NumPol":np.arange(1,len(x),1),"Polinomios": np.transpose(Tabla_string)})
            plt.cla()
            plt.clf()
            for pol, x_el in zip(Tabla, range(len(x)-1)):
                t = np.linspace(x[x_el],x[x_el+1], 100) #domain function
                z=0
                for i, val in zip(range(len(pol)), pol):
                    z += val*t**(len(pol)-i-1)
                plt.plot(t,z)
            plt.title(f"Gráfica de la interpolación de los puntos")
            plt.plot(x,y, marker='.', ls='none', ms=10, color="k")
            plt.legend(["Interpolación","Puntos dados"])
            plt.xlabel("x")
            plt.ylabel("y")
            plt.grid()
            plt.savefig('Grafica.png')
            plt.cla()
            plt.clf()
            context = {'mensaje':mensaje,'df':df,'grafica':'../Grafica.png','x':x,'y':y,'d':d,'paso':True,'tipo':tipo}
            return render(request,template_name='13-spline.html',context=context)

    if request.method == 'GET':
        return render(request,template_name='13-spline.html')