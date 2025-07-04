from django.urls import path
from myapp import views

urlpatterns = [
    path('', views.inicio, name='inicio'),
    path('login/', views.vistaLogin, name='inicio'),
    path('registro/', views.registroUser),
    path('inicio/',views.custom_login),
    path('usAdmin/', views.inicioAdmin, name="dashboard-adm"),
    path('usEstudiante/', views.inicioEstudiante, name="student_dashboard"),
    path('usDocente/', views.inicioDocente, name="teacher_dashboard"),
    path("modificar-molecula/", views.modificar_molecula, name="modificar_molecula"),
    path("editor/", views.vista_modificador, name="vista_modificador"),
    path('editor-dibujo/', views.editor_dibujo, name='editor_dibujo'),
    path('procesar-molecula/', views.procesar_molecula, name='procesar_molecula'),
    path('registroExitoso/',views.registroExitoso),


]