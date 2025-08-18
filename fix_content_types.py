import os
import django

# Configurar Django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mysite.settings')
django.setup()

from django.contrib.contenttypes.models import ContentType
from user.models import PreguntaCuestionario, Cuestionario

def fix_content_types():
    print("🔧 Corrigiendo content_types...")
    
    # Obtener el ContentType de Cuestionario
    cuestionario_ct = ContentType.objects.get_for_model(Cuestionario)
    print(f"✅ ContentType de Cuestionario: {cuestionario_ct.id}")
    
    # Actualizar todas las preguntas existentes
    preguntas_actualizadas = 0
    
    for pregunta in PreguntaCuestionario.objects.all():
        # Si tiene cuestionario (relación vieja)
        if hasattr(pregunta, 'cuestionario') and pregunta.cuestionario:
            pregunta.content_type = cuestionario_ct
            pregunta.object_id = pregunta.cuestionario.id
            pregunta.save()
            preguntas_actualizadas += 1
            print(f"✅ Pregunta {pregunta.id} actualizada para cuestionario {pregunta.cuestionario.id}")
    
    print(f"🎉 {preguntas_actualizadas} preguntas actualizadas exitosamente")

if __name__ == "__main__":
    fix_content_types()
    