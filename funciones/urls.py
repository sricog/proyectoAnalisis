from django.urls import path # Como este es un archivo creado por mi, yo mismo pongo la importación
from funciones import views
from django.conf.urls.static import static
from django.conf import settings


urlpatterns = [
    path('',views.menu_page, name='menu'),
    path('biseccion/',views.biseccion_page, name='biseccion'),
    path('regla_falsa/',views.regla_falsa_page, name='regla_falsa'),
    path('punto_fijo/',views.punto_fijo_page, name='punto_fijo'),
    path('newton/',views.newton_page, name='newton'),
    path('raices_multiples/',views.raices_multiples_page, name='raices_multiples'),
    path('secante/',views.secante_page, name='secante'),
    path('jacobi/',views.jacobi_page, name='jacobi'),
    path('gauss_seidel/',views.gauss_seidel_page, name='gauss_seidel'),
    path('SOR/',views.SOR_page, name='SOR'),
    path('vandermonde/',views.vandermonde_page, name='vandermonde'),
    path('newton_diferencias_divididas/',views.newton_diferencias_divididas_page, name='newton_diferencias_divididas'),
    path('lagrange/',views.lagrange_page, name='lagrange'),
    path('spline/',views.spline_page, name='spline'),
]

urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)