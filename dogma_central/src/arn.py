"""
arn.py
------
Simula la TRANSCRIPCIÓN: ADN -> ARN mensajero (ARNm).

Según el material de la asignatura:
    - La ARN polimerasa (con el factor sigma en bacterias, que reconoce el
      promotor) usa una de las dos hebras de ADN como MOLDE, leyéndola en
      sentido 3' -> 5'.
    - El ARNm se sintetiza en sentido 5' -> 3', complementario a esa hebra
      molde.
    - Se mantiene la complementariedad de bases, sustituyendo la Timina (T)
      por Uracilo (U).
    - El resultado tiene la misma secuencia que la hebra codificante del
      ADN, salvo T -> U.
"""

from __future__ import annotations

from dataclasses import dataclass

from .adn import ADN

# Complementariedad ADN(molde) -> ARN
COMPLEMENTO_ADN_A_ARN = {"A": "U", "T": "A", "C": "G", "G": "C"}


@dataclass
class ResultadoTranscripcion:
    hebra_molde_usada: str       # hebra de ADN 3' <- 5' usada como plantilla
    hebra_codificante: str       # hebra de ADN 5' -> 3' (para comparar con el ARNm)
    arnm: str                    # ARN mensajero resultante, 5' -> 3'


class ARNm:
    """ARN mensajero resultante de transcribir una molécula de ADN."""

    def __init__(self, secuencia: str, nombre: str = "ARNm"):
        self.nombre = nombre
        self.secuencia = secuencia.upper()  # 5' -> 3'

    def __len__(self) -> int:
        return len(self.secuencia)

    def __repr__(self) -> str:
        return f"ARNm('{self.nombre}', {len(self)} nt)"


def transcribir(adn: ADN, nombre: str = "ARNm") -> ResultadoTranscripcion:
    """
    Transcribe la hebra molde de una molécula de ADN a ARN mensajero.

    La hebra molde (`adn.hebra_molde`) está alineada base a base con la
    hebra codificante y se guarda con su extremo 3' a la izquierda y su 5' a
    la derecha. La ARN polimerasa la recorre en sentido 3' -> 5' (es decir,
    de izquierda a derecha tal y como la tenemos guardada) y sintetiza el
    ARNm en 5' -> 3', también de izquierda a derecha.
    """
    molde_3_a_5 = adn.hebra_molde  # ya alineado 3' <- 5' con la codificante
    # Cada base del molde impone, por complementariedad, la base del ARNm
    # que queda exactamente en su misma posición (A-U, T-A, C-G, G-C).
    # El resultado global coincide con la hebra codificante sustituyendo
    # T por U, tal y como indica la teoría de la asignatura.
    arnm_secuencia = "".join(COMPLEMENTO_ADN_A_ARN[base] for base in molde_3_a_5)

    return ResultadoTranscripcion(
        hebra_molde_usada=molde_3_a_5,
        hebra_codificante=adn.hebra_codificante,
        arnm=arnm_secuencia,
    )
