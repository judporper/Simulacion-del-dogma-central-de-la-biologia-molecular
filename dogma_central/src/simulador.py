"""
simulador.py
------------
Clase SimuladorDogmaCentral: une los tres módulos (adn, arn, proteina) y
ejecuta el flujo completo ADN → ADN → ARN → proteína.
"""

from __future__ import annotations

from dataclasses import dataclass

from .adn import ADN, ResultadoReplicacion
from .arn import ResultadoTranscripcion, transcribir
from .proteina import ResultadoTraduccion, traducir


@dataclass
class ResultadoSimulacion:
    adn: ADN
    replicacion: ResultadoReplicacion
    transcripcion: ResultadoTranscripcion
    traduccion: ResultadoTraduccion


class SimuladorDogmaCentral:
    """Orquesta la simulación completa sobre una secuencia de ADN dada."""

    def __init__(self, secuencia_adn: str, nombre: str = "Secuencia sin nombre"):
        self.adn = ADN(secuencia_adn, nombre=nombre)

    def ejecutar(self, tamano_fragmento_okazaki: int = 10) -> ResultadoSimulacion:
        replicacion = self.adn.replicar(tamano_fragmento_okazaki=tamano_fragmento_okazaki)
        transcripcion = transcribir(self.adn)
        traduccion = traducir(transcripcion.arnm)
        return ResultadoSimulacion(
            adn=self.adn,
            replicacion=replicacion,
            transcripcion=transcripcion,
            traduccion=traduccion,
        )
