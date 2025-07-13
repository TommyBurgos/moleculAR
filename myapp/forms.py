from django import forms
from user.models import Curso, Recurso

class CursoForm(forms.ModelForm):
    class Meta:
        model = Curso
        fields = ['titulo', 'descripcion', 'imagen_portada', 'categoria']
        widgets = {
            'titulo': forms.TextInput(attrs={'class': 'form-control'}),
            'descripcion': forms.Textarea(attrs={'class': 'form-control', 'rows': 4}),
            'categoria': forms.Select(attrs={'class': 'form-select'}),
        }

class RecursoForm(forms.ModelForm):
    class Meta:
        model = Recurso
        fields = ['titulo', 'descripcion', 'tipo', 'contenido_texto', 'video_url']
        widgets = {
            'contenido_texto': forms.Textarea(attrs={'rows': 4}),
            'video_url': forms.HiddenInput(),  # Se mostrará con JS al subir el video
        }

    def clean(self):
        cleaned_data = super().clean()
        tipo = cleaned_data.get('tipo')
        contenido = cleaned_data.get('contenido_texto')
        video = cleaned_data.get('video_url')

        if tipo and tipo.nombre.lower() == 'texto' and not contenido:
            raise forms.ValidationError("Debe ingresar el contenido del texto.")
        if tipo and tipo.nombre.lower() == 'video' and not video:
            raise forms.ValidationError("Debe subir un video.")
