from django.shortcuts import render, redirect, get_object_or_404
from django.core.paginator import Paginator
from django.db.models import Q
from django.http import HttpResponse
from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter
from django.middleware.csrf import rotate_token
from django.contrib.auth.decorators import login_required
from django.contrib import messages
import tempfile
import json


from user.models import User,Rol, Curso, Seccion, TipoRecurso, Recurso, Cuestionario, Practica, Modelo, InscripcionCurso, MoleculaEstudiante, Competencia, PreguntaCuestionario,TipoPregunta, ProgresoUsuario
from django.contrib.auth import login, logout, authenticate
import re
from .permissions import role_required

@csrf_exempt
def modificar_molecula(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        smiles = data.get("smiles")
        atomo_idx = int(data.get("atomo_idx"))
        grupo = data.get("grupo")

        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        edmol = Chem.RWMol(mol)

        try:
            # Remover un hidrógeno unido al átomo objetivo
            h_idx = None
            for vecino in edmol.GetAtomWithIdx(atomo_idx).GetNeighbors():
                if vecino.GetAtomicNum() == 1:
                    h_idx = vecino.GetIdx()
                    break

            if h_idx is not None:
                edmol.RemoveAtom(h_idx)

            if grupo == "CH3":
                c = Chem.Atom(6)
                h1 = Chem.Atom(1)
                h2 = Chem.Atom(1)
                h3 = Chem.Atom(1)
                c_idx = edmol.AddAtom(c)
                edmol.AddBond(atomo_idx, c_idx, Chem.rdchem.BondType.SINGLE)
                edmol.AddBond(c_idx, edmol.AddAtom(h1), Chem.rdchem.BondType.SINGLE)
                edmol.AddBond(c_idx, edmol.AddAtom(h2), Chem.rdchem.BondType.SINGLE)
                edmol.AddBond(c_idx, edmol.AddAtom(h3), Chem.rdchem.BondType.SINGLE)

            mol_final = edmol.GetMol()
            Chem.SanitizeMol(mol_final)
            mol_final.UpdatePropertyCache(strict=False)
            AllChem.EmbedMolecule(mol_final)
            AllChem.UFFOptimizeMolecule(mol_final)

            # Guardar SDF temporal
            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp:
                writer = SDWriter(temp.name)
                writer.write(mol_final)
                writer.close()
                with open(temp.name, "r") as f:
                    sdf_data = f.read()

            return JsonResponse({"sdf": sdf_data})
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=400)

    # ✅ Manejo para GET
    return JsonResponse({"message": "Usa POST para modificar la molécula."})

def vista_modificador(request):
    return render(request, 'modificador_molecula.html')

# Create your views here.
def inicio(request):
    return render(request, 'index.html')

def registroUser(request):
    return render(request, 'sign-up.html')

#REGISTRAR USUARIO ESTUDIANTE
def registroExitoso(request):
    if request.method == 'POST':
        print(request.POST)
        if request.POST['password'] == request.POST['confirmPassword']:
            try:
                password = request.POST.get('password')
                email = request.POST.get('email')                
                # Validaciones
                if len(password) < 8 or not re.match(r'^(?=.*[A-Za-z])(?=.*\d)[A-Za-z\d]{8,}$', password):
                    return HttpResponse('La contraseña debe tener al menos 8 caracteres y contener letras y números.')

                if not re.match(r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$', email):
                    return HttpResponse('El formato del correo electrónico no es válido.')
                
                # Crear usuario
                user = User.objects.create_user(
                    username=email,  # Usa el correo como username
                    email=email,
                    password=password,
                    first_name=' ',
                    last_name=' '
                )
                print('***INICIO***')
                # Si quieres asignar rol por defecto
                rol_estudiante = Rol.objects.get(nombre='Estudiante')
                print(rol_estudiante)
                user.rol = rol_estudiante
                print(user.rol)                                
                user.save()
                print('***FIN***')
                login(request, user)
                print("Usuario creado y logueado exitosamente.")
                
                return redirect('inicio')  # O la vista que tengas definida
            except Exception as e:
                print(f"❌ Error al crear el usuario: {e}")
                return HttpResponse('Ocurrió un error al crear el usuario.')
        else:
            return HttpResponse('Las contraseñas no coinciden.')
    return HttpResponse('Método no permitido')


def vistaLogin(request):
    return render(request,'sign-in.html')

@role_required('Admin')
def inicioAdmin(request):
    user = request.user    
    imgPerfil=user.imgPerfil  
    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario':user.username,        
    }   
    return render(request, 'usAdmin/admin-dashboard.html', context)

@role_required('Admin')
def vistaAllCursos(request):
    user = request.user    
    imgPerfil=user.imgPerfil 
    cursos = Curso.objects.all().order_by('-fecha_creacion') 
    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario':user.username, 
        'cursos': cursos       
    }   
    return render(request,'usAdmin/detalleCursos.html',context)
from django.db.models import Prefetch

@role_required('Admin')
@login_required
def detalle_curso(request, curso_id):
    curso = get_object_or_404(Curso, id=curso_id)
    
    recursos_prefetch = Prefetch(
        'recursos',
        queryset=Recurso.objects.order_by('orden'),
        to_attr='recursos_ordenados'
    )
    
    secciones = Seccion.objects.filter(curso=curso_id).prefetch_related(recursos_prefetch).order_by('orden')

    recursosTipo = TipoRecurso.objects.all()
    user = request.user
    imgPerfil = user.imgPerfil

    return render(request, 'usAdmin/detalle_curso.html', {
        'curso': curso,
        'imgPerfil': imgPerfil,
        'usuario': user.username,
        'secciones': secciones,
        'tipos_recurso': recursosTipo
    })

from django.utils.html import strip_tags

def limpiar_html_permitiendo_alineacion(html):
    if not html:
        return ""

    # Elimina todo menos etiquetas <p>, <strong>, <em>, <u>, <br>, <span>
    html = re.sub(r'<(?!/?(p|strong|em|u|br|span|img)(\s|>|/))[^>]+>', '', html)

    # Filtra los estilos permitidos (solo text-align)
    html = re.sub(r'style="[^"]*"', lambda m: filtrar_estilos(m.group()), html)

    return html

def filtrar_estilos(style_str):
    estilos_permitidos = ['text-align', 'color', 'background-color']
    estilos = style_str.replace('style="', '').replace('"', '').split(';')
    estilos_filtrados = [s for s in estilos if any(permitido in s for permitido in estilos_permitidos)]
    if estilos_filtrados:
        return f'style="{";".join(estilos_filtrados)}"'
    return ''


import bleach

@login_required
def detalle_recurso(request, recurso_id):
    print("He ingresado al detalle del recurso")
    recurso = get_object_or_404(Recurso, id=recurso_id)

    if request.user.is_authenticated and hasattr(request.user, 'rol') and request.user.rol.nombre == 'Estudiante':
        progreso, creado = ProgresoUsuario.objects.get_or_create(
            usuario=request.user,
            recurso=recurso,
            defaults={'visto': True}
        )
        if not creado and not progreso.visto:
            progreso.visto = True
            progreso.save()

    curso = recurso.seccion.curso 
    contenido_html = recurso.contenido_texto or ""
    html_limpio = limpiar_html_permitiendo_alineacion(contenido_html)

    recursos_prefetch = Prefetch(
        'recursos',
        queryset=Recurso.objects.order_by('orden'),
        to_attr='recursos_ordenados'
    )    
    secciones = Seccion.objects.filter(curso=curso).prefetch_related(recursos_prefetch).order_by('orden')

    return render(request, 'usAdmin/detalleRecurso.html', {
        'recurso': recurso,
        'secciones': secciones,
        'curso': curso,
        'html_limpio': html_limpio
    })

@role_required('Admin')
@login_required
def editar_video(request, recurso_id):
    print("Acabo de ingresar a editar texto")
    recurso = get_object_or_404(Recurso, id=recurso_id, tipo__nombre='Video')

    if request.method == 'POST':
        recurso.titulo = request.POST.get('titulo')
        recurso.descripcion = request.POST.get('descripcion')
        # Opcionalmente permitir reemplazar el video
        if 'archivo' in request.FILES:
            recurso.archivo = request.FILES['archivo']
        recurso.save()
        messages.success(request, "Video actualizado correctamente.")
        return redirect('detalle_curso', curso_id=recurso.seccion.curso.id)

    return render(request, 'usAdmin/editar_video.html', {'recurso': recurso})


@role_required('Admin')
@login_required
def editar_texto(request, recurso_id):
    print("Acabo de ingresar a editar texto")
    recurso = get_object_or_404(Recurso, id=recurso_id, tipo__nombre='Texto')
    print(recurso)
    if request.method == 'POST':
        recurso.titulo = request.POST.get('titulo')
        recurso.descripcion = request.POST.get('descripcion')
        recurso.contenido_texto = request.POST.get('contenido_texto')
        recurso.save()
        messages.success(request, "Texto actualizado correctamente.")
        return redirect('detalle_curso', curso_id=recurso.seccion.curso.id)

    return render(request, 'usAdmin/editar_texto.html', {'recurso': recurso})



def crear_modulo(request, curso_id):
    curso = get_object_or_404(Curso, id=curso_id)

    if request.method == 'POST':
        titulo = request.POST.get('titulo', '').strip()
        descripcion = request.POST.get('descripcion', '').strip()
        orden = request.POST.get('orden', 1)

        if not titulo:
            messages.error(request, "El título del módulo es obligatorio.")
            return redirect('crear_modulo', curso_id=curso.id)

        Seccion.objects.create(
            curso=curso,
            titulo=titulo,
            descripcion=descripcion,
            orden=orden
        )
        messages.success(request, "Módulo creado exitosamente.")
        return redirect('detalle_curso', curso_id=curso.id)

    return render(request, 'usAdmin/detalleCursos.html', {'curso': curso})

#FUNCIÓN QUE ME AYUDA A OBTENER EL SMILE
import requests
def obtener_smiles_pubchem(nombre_molecula):
    nombre_molecula = nombre_molecula.strip().lower()
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{nombre_molecula}/property/IsomericSMILES,SMILES/JSON"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            properties = data.get('PropertyTable', {}).get('Properties', [])
            if properties:
                return (
                    properties[0].get('IsomericSMILES') or
                    properties[0].get('SMILES')
                )
    except Exception as e:
        print(f"❌ Error al consultar PubChem: {e}")
    return None

def ver_practica(request, recurso_id):
    recurso = get_object_or_404(Recurso, id=recurso_id)
    practica = get_object_or_404(Practica, recurso=recurso)
    print(f"el recurso es {recurso}")
    print(f"la practica es {practica}")

    modelo = practica.modelo_objetivo
    print(f"el modelo es {modelo}")
    print(f"el smile es {modelo.smiles}")
    print(modelo.smiles if modelo else "")
    smiles = modelo.smiles if modelo else ""
    
    return render(request, 'usAdmin/ver_practica.html', {
        'practica': practica,
        'recurso': recurso,
        'modelo':modelo,
        'smiles': smiles
    })


@login_required
def editar_practica(request, recurso_id):
    recurso = get_object_or_404(Recurso, id=recurso_id)
    practica, creado = Practica.objects.get_or_create(recurso=recurso)

    if request.method == 'POST':
        fuente = request.POST.get('fuente_molecula')
        titulo_objetivo = request.POST.get('titulo_objetivo')
        instrucciones = request.POST.get('instrucciones')
        smiles = None
        modelo = None

        if fuente == 'jsme':
            smiles = request.POST.get('mol_input')
            if smiles and smiles.strip():
                modelo = Modelo.objects.create(
                    titulo=titulo_objetivo,
                    descripcion=instrucciones,
                    tipo='2D',
                    smiles=smiles.strip(),
                    creado_por=request.user,
                    visible_biblioteca=False
                )
            else:
                return render(request, 'usAdmin/editor_dibujo.html', {
                    'recurso': recurso,
                    'error': 'No se pudo obtener el SMILES desde el editor JSME.'
                })

        elif fuente == 'pubchem':
            nombre_molecula = request.POST.get('nombre_molecula')
            print("Nombre recibido desde formulario:", nombre_molecula)
            smiles = obtener_smiles_pubchem(nombre_molecula)
            print("SMILES recibido desde PubChem:", smiles)
            if smiles and smiles.strip():
                modelo = Modelo.objects.create(
                    titulo=nombre_molecula,
                    descripcion='Molécula importada desde PubChem',
                    tipo='2D',
                    smiles=smiles.strip(),
                    creado_por=request.user,
                    visible_biblioteca=False
                )
            else:
                return render(request, 'usAdmin/editor_dibujo.html', {
                    'recurso': recurso,
                    'error': f"No se pudo encontrar un SMILES válido para '{nombre_molecula}'. Intenta con otro nombre."
                })

        # Crear o actualizar la práctica solo si se creó un modelo
        if modelo:
            practica, creada = Practica.objects.get_or_create(
                recurso=recurso,
                defaults={
                    'titulo_objetivo': titulo_objetivo,
                    'instrucciones': instrucciones,
                    'modelo_objetivo': modelo
                }
            )

            if not creada:
                practica.titulo_objetivo = titulo_objetivo
                practica.instrucciones = instrucciones
                practica.modelo_objetivo = modelo
                practica.save()

            return redirect('detalle_curso', curso_id=recurso.seccion.curso.id)
        else:
            return render(request, 'usAdmin/editor_dibujo.html', {
                'recurso': recurso,
                'error': 'No se pudo crear la práctica porque no se generó ningún modelo.'
            })

    # Enviar el SMILES existente si hay modelo asociado, si no, un string vacío
    smiles_existente = practica.modelo_objetivo.smiles if practica.modelo_objetivo else ""
    print(f"el smile existente es {smiles_existente}")

    return render(request, 'usAdmin/editor_dibujo.html', {
        'recurso': recurso,
        'practica': practica,
        'instruccion': practica.instrucciones if practica.instrucciones else "",
        'smiles': smiles_existente
    })


def crear_recurso(request, seccion_id):
    if request.method == 'POST':
        print("Entré al post de crear recurso")
        # Datos básicos del recurso
        titulo = request.POST.get('titulo')
        print(f"Titulo: {titulo}")
        descripcion = request.POST.get('descripcion')
        imagen = request.FILES.get('imagen_referencial')
        tipo_id = request.POST.get('tipo')
        visible = 'visible_biblioteca' in request.POST
        print("es visible? ", visible)
        seccion_id = request.POST.get('seccion_id')
        print("La sección id es" + seccion_id)
        tipo = TipoRecurso.objects.get(id=tipo_id)
        seccion = Seccion.objects.get(id=seccion_id)
        contenido_texto = request.POST.get('contenido_texto')

        recurso = Recurso.objects.create(
            titulo=titulo,
            descripcion=descripcion,
            tipo=tipo,
            imagen=imagen,
            seccion=seccion,
            visible_biblioteca=visible,
            contenido_texto=contenido_texto
        )
        print(f"Recurso creado {recurso}")
        # CREACIÓN SEGÚN TIPO DE RECURSO
        nombre_tipo = tipo.nombre.lower()
        print(f"nombre del tipo de recurso es {nombre_tipo}")
        print('video' in nombre_tipo)
        print("Comparando si ingresa a pràctica")
        print('Practica' in nombre_tipo)
        print('practica' in nombre_tipo)
        print(tipo.nombre.lower() == 'practica')

        if 'video' in nombre_tipo:
            print("Dentro del if")
            video_url = request.POST.get('video_url')
            print(video_url)
            if video_url:
                print(recurso.video_url)
                recurso.video_url = video_url
                print(recurso.video_url)
                recurso.save()
            else:
                print("⚠️ No se recibió la URL del video desde el formulario.")


        elif 'cuestionario' in nombre_tipo:
            print("Dentro del if cuestionario")
            instrucciones = request.POST.get('nombre_cuestionario')
            tiempo_limite = request.POST.get('descripcion_cuestionario')
            Cuestionario.objects.create(
                recurso=recurso,
                instrucciones=instrucciones,
                tiempo_limite=tiempo_limite
            )        
        elif tipo.nombre.lower() == 'practica':
            # Por ahora NO creamos la práctica. Solo dejamos el recurso creado.
            # Luego redirigiremos a una vista donde se edita ese recurso.
            print("Se creo el recurso practica")
            return redirect('editar_practica', recurso.id)  
        
        elif tipo.nombre.lower() == 'competencia':
            print("Dentro del if competencia")
            instrucciones = request.POST.get('instrucciones_competencia', '')
            # Aquí puedes generar un código aleatorio o permitir ingresarlo desde el formulario
            import random, string
            codigo = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))

            Competencia.objects.create(
                recurso=recurso,
                instrucciones=instrucciones,
                codigo_acceso=codigo
            )
            print("Competencia creada correctamente.")
            return redirect('panel_competencia', recurso.competencia.id)

        return redirect('detalle_recurso',recurso.id)

@role_required('Admin')
def detalleUsuarios(request):
    user = request.user
    usuarios2 = User.objects.all().order_by('-date_joined')
    usuarios = User.objects.all().order_by('-date_joined')
    imgPerfil=user.imgPerfil

    usuarios_con_progreso = []    
    busqueda = request.GET.get("buscar")
    if busqueda:
        usuarios = User.objects.filter(
            Q(username__icontains=busqueda) |
            Q(first_name__icontains=busqueda) |
            Q(last_name__icontains=busqueda) |
            Q(email__icontains=busqueda)
        ).distinct()
        print("entre al if de busqueda")
    print("LISTADO DE USUARIOS")
    print(usuarios)
    print(usuarios2)

        # Paginación — debe ir después de filtrar los usuarios
    paginator = Paginator(usuarios, 10)  # 10 usuarios por página
    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)
    # 1. Cantidad de usuarios cuyo rol es igual a 2  
    
    filtro_progreso = request.GET.get('progreso')

    for usuario in page_obj:
        progreso = 0
        recursos_totales = 0
        recursos_vistos = 0
        print("viendo viendoo")
        print(hasattr(usuario, 'rol') and usuario.rol == 'Estudiante')
        print(usuario.rol)        
        print(str(usuario.rol) == 'Estudiante')
        print(hasattr(usuario, 'rol'))

        if hasattr(usuario, 'rol') and str(usuario.rol) == 'Estudiante':
            # Cursos del estudiante
            cursos = Curso.objects.filter(inscripcioncurso__estudiante=usuario)
            recursos_totales = Recurso.objects.filter(seccion__curso__in=cursos).count()
            recursos_vistos = ProgresoUsuario.objects.filter(usuario=usuario, visto=True).count()
            progreso = int((recursos_vistos / recursos_totales) * 100) if recursos_totales > 0 else 0
        # Aplicamos el filtro si se solicitó
        if filtro_progreso == 'bajo' and progreso >= 50:
            continue  # ignorar este usuario
        usuarios_con_progreso.append({
            'usuario': usuario,
            'progreso': progreso,
            'recursos_vistos': recursos_vistos,
            'recursos_totales': recursos_totales,
        }) 
      
    context = {                     
        'imgPerfil': imgPerfil,        
        'usuario':user.username,  
        'usuarios':page_obj,
        'usuarios_con_progreso': usuarios_con_progreso,  # Lista con progreso incluido
      
    }
    return render(request, 'usAdmin/detalle_usuarios.html', context)

@role_required('Admin')
def editar_usuario(request, user_id):
    usuario = get_object_or_404(User, id=user_id)
    cursos_inscritos = InscripcionCurso.objects.filter(estudiante=usuario).select_related('curso')
    # Total de recursos en la plataforma o en los cursos donde está inscrito
    cursos = Curso.objects.filter(inscripcioncurso__estudiante=usuario)
    recursos_totales = Recurso.objects.filter(seccion__curso__in=cursos).count()

    # Recursos vistos por el estudiante
    recursos_vistos = ProgresoUsuario.objects.filter(usuario=usuario, visto=True).count()

    # Evitamos división por cero
    if recursos_totales > 0:
        progreso = int((recursos_vistos / recursos_totales) * 100)
    else:
        progreso = 0

    if request.method == 'POST':
        usuario.first_name = request.POST.get('first_name')
        usuario.last_name = request.POST.get('last_name')
        usuario.email = request.POST.get('email')
        usuario.is_active = 'is_active' in request.POST  # checkbox activo
        usuario.save()
        messages.success(request, 'Usuario actualizado correctamente.')
        return redirect('detalleUsuarios-adm')  # o donde quieras redirigir

    context = {
        'usuario': usuario,
        'cursos_inscritos': cursos_inscritos,
        'progreso': progreso,
        'recursos_totales': recursos_totales,
        'recursos_vistos': recursos_vistos,

    }
    return render(request, 'usAdmin/editar_usuario.html', context)

@role_required('Admin')
def eliminar_usuario(request, user_id):
    usuario = get_object_or_404(User, id=user_id)
    if request.method == 'POST':
        usuario.delete()
        messages.success(request, "Usuario eliminado correctamente.")
        return redirect('detalle_usuarios')
    return render(request, 'usAdmin/confirmar_eliminar_usuario.html', {'usuario': usuario})


@role_required('Admin')
@login_required
def resetear_contrasena(request, id):
    usuario = get_object_or_404(User, id=id)
    
    return redirect('detalle_usuarios')

@role_required('Admin')
def resetContra_usuario(request, id):
    usuario = get_object_or_404(User, id=id)
    # lógica para mostrar o editar
    nueva_contrasena = "12345678"
    usuario.set_password(nueva_contrasena)  # Cambia la contraseña de forma segura
    usuario.save()
    messages.success(request, f"La contraseña del usuario '{usuario.username}' ha sido reseteada con éxito.")
    return redirect('detalleUsuarios-adm')




import boto3
import uuid
from django.conf import settings

@csrf_exempt
def obtener_presigned_url(request):
    if request.method == 'POST':
        nombre_archivo = request.POST.get('filename')
        extension = nombre_archivo.split('.')[-1]
        nuevo_nombre = f"videos/{uuid.uuid4()}.{extension}"

        s3 = boto3.client(
            's3',
            aws_access_key_id=settings.AWS_ACCESS_KEY_ID,
            aws_secret_access_key=settings.AWS_SECRET_ACCESS_KEY,
            region_name=settings.AWS_S3_REGION_NAME
        )

        # ⚠️ NO incluir campos ACL ni condiciones ACL
        presigned_post = s3.generate_presigned_post(
            Bucket=settings.AWS_STORAGE_BUCKET_NAME,
            Key=nuevo_nombre,
            Fields={},  # sin ACL aquí
            Conditions=[],
            ExpiresIn=3600
        )

        url_publica = f"https://{settings.AWS_STORAGE_BUCKET_NAME}.s3.{settings.AWS_S3_REGION_NAME}.amazonaws.com/{nuevo_nombre}"

        return JsonResponse({
            'data': presigned_post,
            'url': url_publica
        })





from .forms import CursoForm

@role_required('Admin')
def vistaCrearCurso(request):
    user = request.user    
    imgPerfil=user.imgPerfil 
    if request.method == 'POST':
        form = CursoForm(request.POST, request.FILES)
        if form.is_valid():
            curso = form.save(commit=False)
            curso.profesor = request.user  # Asignar el usuario actual como profesor
            curso.save()
            return redirect('detalleCursos-adm', curso.id) 
    else:
        form = CursoForm() 
    context = {                  
        'imgPerfil': imgPerfil,        
        'usuario':user.username,
        'form': form        
    }   
    return render(request,'usAdmin/crearCurso.html',context)

#EMPEZANDO INTENTAR CREAR COMPETENCIAS

def panel_competencia(request, competencia_id):
    competencia = get_object_or_404(Competencia, id=competencia_id)
    preguntas = competencia.preguntas.order_by('orden')

    return render(request, 'usAdmin/panel_control.html', {
        'competencia': competencia,
        'preguntas': preguntas,
    })

def agregar_pregunta_competencia(request, competencia_id):
    competencia = Competencia.objects.get(id=competencia_id)

    if request.method == 'POST':
        enunciado = request.POST.get('enunciado')
        tipo_id = request.POST.get('tipo_pregunta')
        orden = request.POST.get('orden')
        puntaje = request.POST.get('puntaje')

        tipo = TipoPregunta.objects.get(id=tipo_id)
        tipo_nombre = tipo.nombre.lower()

        # Crear la pregunta base
        pregunta = PreguntaCuestionario.objects.create(
            competencia=competencia,
            enunciado=enunciado,
            tipo=tipo,
            orden=orden,
            puntaje=puntaje
        )

        # Opciones según el tipo
        if tipo_nombre == 'seleccion_multiple':
            for i in range(1, 5):
                texto = request.POST.get(f'opcion_{i}')
                es_correcta = request.POST.get('correcta') == str(i)
                if texto:
                    pregunta.opciones.create(texto=texto, es_correcta=es_correcta)

        elif tipo_nombre == 'verdadero_falso':
            respuesta = request.POST.get('vf_respuesta')
            if respuesta:
                pregunta.opciones.create(texto='Verdadero', es_correcta=(respuesta == 'verdadero'))
                pregunta.opciones.create(texto='Falso', es_correcta=(respuesta == 'falso'))

        elif tipo_nombre == 'completar':
            palabra_correcta = request.POST.get('palabra_correcta')
            if palabra_correcta:
                pregunta.opciones.create(texto=palabra_correcta.strip(), es_correcta=True)

        return redirect('panel_competencia', competencia.id)

    tipos = TipoPregunta.objects.all()
    return render(request, 'usAdmin/agregar_pregunta.html', {
        'competencia': competencia,
        'tipos': tipos
    })


from django.http import HttpResponseRedirect
from django.urls import reverse

def iniciar_competencia(request, competencia_id):
    competencia = get_object_or_404(Competencia, id=competencia_id)
    competencia.activa = True
    competencia.save()
    return HttpResponseRedirect(reverse('panel_competencia', args=[competencia.id]))

def enviar_pregunta(request, pregunta_id):
    pregunta = get_object_or_404(PreguntaCuestionario, id=pregunta_id)
    competencia = pregunta.competencia
    competencia.pregunta_actual = pregunta
    competencia.save()
    return HttpResponseRedirect(reverse('panel_competencia', args=[competencia.id]))

##HASTA AQU[I]

@role_required('Estudiante')
def inicioEstudiante(request):
    usuario = request.user
    imgPerfil=usuario.imgPerfil
    inscripciones = InscripcionCurso.objects.filter(estudiante=usuario).select_related('curso')
    cursos = [ins.curso for ins in inscripciones]
    
    print(imgPerfil) 
    context = {       
        'cursos': cursos,           
        'imgPerfil': imgPerfil,        
        'usuario':usuario,                
    }      
    return render(request, 'estudiante/student-dashboard.html', context)    

from django.contrib import messages

def unirse_curso(request):
    if request.method == 'POST':
        codigo = request.POST.get('codigo_acceso')
        try:
            curso = Curso.objects.get(codigo_acceso=codigo)
            ya_inscrito = InscripcionCurso.objects.filter(curso=curso, estudiante=request.user).exists()
            if ya_inscrito:
                messages.warning(request, 'Ya estás inscrito en este curso.')
            else:
                InscripcionCurso.objects.create(curso=curso, estudiante=request.user)
                messages.success(request, 'Te has unido correctamente al curso.')
        except Curso.DoesNotExist:
            messages.error(request, 'El código ingresado no es válido.')
        return redirect('listar_cursos')  # o tu página de cursos

def listar_cursos(request):
    usuario = request.user
    inscripciones = InscripcionCurso.objects.filter(estudiante=usuario).select_related('curso')
    cursos = [ins.curso for ins in inscripciones]

    return render(request, 'estudiante/student-course-list.html', {
        'cursos': cursos,
        'usuario': usuario,
        'imgPerfil': usuario.imgPerfil
    })



def detalle_cursoEstudiante(request, curso_id):
    print("Ingrese al detalle del curso estudiante")
    curso = get_object_or_404(Curso, id=curso_id)
    
    recursos_prefetch = Prefetch(
        'recursos',
        queryset=Recurso.objects.order_by('orden'),
        to_attr='recursos_ordenados'
    )
    print(recursos_prefetch)
    secciones = Seccion.objects.filter(curso=curso_id).prefetch_related(recursos_prefetch).order_by('orden')    
    recursosTipo = TipoRecurso.objects.all()
    user = request.user
    imgPerfil = user.imgPerfil

    return render(request, 'estudiante/student-course-resume.html', {
        'curso': curso,
        'imgPerfil': imgPerfil,
        'usuario': user.username,
        'secciones': secciones,
        'tipos_recurso': recursosTipo
    })

@login_required
def resolver_practica(request, practica_id):
    practica = get_object_or_404(Practica, id=practica_id)

    if request.method == 'POST':
        smiles_estudiante = request.POST.get('mol_input')
        modelo_objetivo = practica.modelo_objetivo

        # Comprobación (simple, se puede mejorar con RDKit)
        es_correcta = smiles_estudiante.strip() == modelo_objetivo.smiles.strip()

        MoleculaEstudiante.objects.create(
            estudiante=request.user,
            practica=practica,
            smiles=smiles_estudiante,
            es_correcta=es_correcta
        )

        return render(request, 'estudiante/resultado_practica.html', {
            'practica': practica,
            'es_correcta': es_correcta
        })

    return render(request, 'estudiante/resolver_practica.html', {
        'practica': practica,
    })

@login_required
def historial_intentos(request, practica_id):
    practica = get_object_or_404(Practica, id=practica_id)
    intentos = MoleculaEstudiante.objects.filter(
        estudiante=request.user,
        practica=practica
    ).order_by('-fecha_envio')

    return render(request, 'estudiante/historial_intentos.html', {
        'practica': practica,
        'intentos': intentos
    })



def biblioteca_practicas(request):
    query = request.GET.get('q', '')  # Capturamos el término de búsqueda

    practicas_visibles = Recurso.objects.filter(
        tipo__nombre__iexact='practica',
        visible_biblioteca=True
    )

    if query:
        practicas_visibles = practicas_visibles.filter(
            Q(titulo__icontains=query)
        )

    # Paginación: 10 por página
    paginator = Paginator(practicas_visibles.order_by('-fecha_creacion'), 10)
    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)

    context = {
        'page_obj': page_obj,
        'query': query
    }
    return render(request, 'estudiante/biblioteca_practicas.html', context)


def inicioDocente(request):
    return render(request, 'docente/instructor-dashboard.html')


def custom_login(request):
    print("entre a la funcion custom")
    if request.method == "POST":  
        print("Estoy en post")   
        print(request.POST) 
        print(request.POST['email'], request.POST['password'])  
        user = authenticate(request, username=request.POST['email'], password= request.POST
        ['password'])
        print(user)
        print("Ya autentique el usuario")
        print(user is not None)
        if user is not None:            
            login(request, user)
            print("Loguie al usuario")
            rol = user.rol_id  # Asignar rol después de la autenticación   
            print(rol)         
            if rol == 1:
                print("Ingrese al rol 1")
                return redirect('dashboard-adm')
            elif rol == 2:
                return redirect('student_dashboard')
            elif rol == 3:
                return redirect('teacher_dashboard')
        else:
            # Manejar error de autenticación
            return render(request, 'generales/accesoDenegado.html', {'error': 'Credenciales inválidas'})
    return redirect('inicio')

def signout(request):
    logout(request)
    rotate_token(request)  # Gira el token CSRF para la nueva sesión
    return redirect('inicio')

from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
import tempfile

@csrf_exempt
def editor_dibujo(request):
    return render(request, 'editor_dibujo.html')

@csrf_exempt
def procesar_molecula(request):
    if request.method == 'POST':
        entrada = request.POST.get('mol_input')
        tipo = request.POST.get('tipo')  # "smiles" o "mol"

        try:
            if tipo == 'smiles':
                mol = Chem.MolFromSmiles(entrada)
                mol = Chem.AddHs(mol)
            elif tipo == 'mol':
                mol = Chem.MolFromMolBlock(entrada)
            else:
                return JsonResponse({'error': 'Formato no válido'}, status=400)

            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)

            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp:
                writer = SDWriter(temp.name)
                writer.write(mol)
                writer.close()
                with open(temp.name, "r") as f:
                    sdf_data = f.read()

            return render(request, 'visor_3d.html', {'sdf_data': sdf_data})

        except Exception as e:
            return JsonResponse({'error': str(e)}, status=400)

    return JsonResponse({'message': 'Usa POST con SMILES o MOL'})
