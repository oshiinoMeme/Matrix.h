#include "s21_matrix.h"

int s21_check_matrix(matrix_t *A) {
  int res = 1;
  if (!A || (*A).rows <= 0 || (*A).columns <= 0 || abs((*A).rows) == INFINITY ||
      abs((*A).columns) == INFINITY || (*A).rows == NAN || (*A).columns == NAN)
    res = 0;
  return res;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = INCORRECT_MATRIX;
  if (rows > 0 && columns > 0) {
    error = OK;
    (*result).rows = rows;
    (*result).columns = columns;
    (*result).matrix = (double **)malloc(sizeof(double *) * rows);
    if ((*result).matrix != NULL) {
      for (int i = 0; i < rows; i++) {
        (*result).matrix[i] = (double *)malloc(sizeof(double) * columns);
      }
      for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
          (*result).matrix[i][j] = 0;
        }
      }
    } else {
      error = INCORRECT_MATRIX;
    }
  } else {
    error = INCORRECT_MATRIX;
  }
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if ((*A).matrix != NULL) {
    for (int i = 0; i < (*A).rows; i++) {
      if ((*A).matrix[i] != NULL) {
        free((*A).matrix[i]);
      }
    }
    free((*A).matrix);
    (*A).matrix = NULL;
  }
  if ((*A).rows) {
    (*A).rows = 0;
  }
  if ((*A).columns) {
    (*A).columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int result = FAILURE;
  if (s21_check_matrix(A) && s21_check_matrix(B)) {
    if ((*A).rows == (*B).rows && (*A).columns == (*B).columns) {
      result = SUCCESS;
      for (int i = 0; i < (*A).rows; i++) {
        for (int j = 0; j < (*A).columns; j++) {
          if (fabs((*A).matrix[i][j] - (*B).matrix[i][j]) >= S21_EPS) {
            result = FAILURE;
          }
        }
      }
    }
  }
  return result;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = INCORRECT_MATRIX;
  if (s21_check_matrix(A) && s21_check_matrix(B) && result != NULL) {
    error = CALC_ERROR;
    if ((*A).rows == (*B).rows && (*A).columns == (*B).columns) {
      error = s21_create_matrix((*A).rows, (*A).columns, result);
      for (int i = 0; i < (*A).rows; i++) {
        for (int j = 0; j < (*A).columns; j++) {
          (*result).matrix[i][j] = (*A).matrix[i][j] + (*B).matrix[i][j];
        }
      }
    }
  }
  return error;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = INCORRECT_MATRIX;
  if (s21_check_matrix(A) && s21_check_matrix(B) && result != NULL) {
    error = CALC_ERROR;
    if ((*A).rows == (*B).rows && (*A).columns == (*B).columns) {
      error = s21_create_matrix((*A).rows, (*A).columns, result);
      for (int i = 0; i < (*A).rows; i++) {
        for (int j = 0; j < (*A).columns; j++) {
          (*result).matrix[i][j] = (*A).matrix[i][j] - (*B).matrix[i][j];
        }
      }
    }
  }
  return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = INCORRECT_MATRIX;
  if (s21_check_matrix(A)) {
    error = s21_create_matrix((*A).rows, (*A).columns, result);
    for (int i = 0; i < (*A).rows; i++) {
      for (int j = 0; j < (*A).columns; j++) {
        (*result).matrix[i][j] = (*A).matrix[i][j] * number;
      }
    }
  }
  return error;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = INCORRECT_MATRIX;
  if (s21_check_matrix(A) && s21_check_matrix(B)) {
    error = CALC_ERROR;
    if ((*A).columns == (*B).rows) {
      error = s21_create_matrix((*A).rows, (*B).columns, result);
      for (int i = 0; i < (*A).rows; i++) {
        for (int j = 0; j < (*B).columns; j++) {
          for (int k = 0; k < (*B).rows; k++) {
            (*result).matrix[i][j] += (*A).matrix[i][k] * (*B).matrix[k][j];
          }
        }
      }
    }
  }
  return error;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error = INCORRECT_MATRIX;
  if (s21_check_matrix(A)) {
    error = s21_create_matrix((*A).columns, (*A).rows, result);
    for (int i = 0; i < (*A).columns; i++) {
      for (int j = 0; j < (*A).rows; j++) {
        (*result).matrix[i][j] = (*A).matrix[j][i];
      }
    }
  }
  return error;
}

void s21_minor_for_det(matrix_t *matrix, int row, int col, matrix_t *minor) {
  int row_minor = 0;
  for (int i = 0; i < row - 1; i++) {
    if (i == 0) {
      row_minor = 1;
    }
    int col_minor = 0;
    for (int j = 0; j < row - 1; j++) {
      if (j == col) {
        col_minor = 1;
      }
      (*minor).matrix[i][j] = (*matrix).matrix[i + row_minor][j + col_minor];
    }
  }
}

double s21_det_calc(matrix_t *matrix) {
  double det = 0;
  int degree = 1;
  if ((*matrix).rows == 1) {
    det = (*matrix).matrix[0][0];
  } else if ((*matrix).rows == 2) {
    det = (*matrix).matrix[0][0] * (*matrix).matrix[1][1] -
          (*matrix).matrix[0][1] * (*matrix).matrix[1][0];
  } else {
    matrix_t minor_matr;
    s21_create_matrix((*matrix).rows - 1, (*matrix).rows - 1, &minor_matr);
    for (int j = 0; j < (*matrix).rows; j++) {
      s21_minor_for_det(matrix, (*matrix).rows, j, &minor_matr);
      det = det + (degree * (*matrix).matrix[0][j] * s21_det_calc(&minor_matr));
      degree = -degree;
    }
    s21_remove_matrix(&minor_matr);
  }
  return det;
}

int s21_determinant(matrix_t *A, double *result) {
  int error = INCORRECT_MATRIX;
  if (s21_check_matrix(A)) {
    error = CALC_ERROR;
    if ((*A).rows == (*A).columns) {
      error = OK;
      *result = (*A).matrix[0][0];
      if ((*A).rows != 1) {
        *result = s21_det_calc(A);
      }
    }
  }
  return error;
}

double s21_det_for_minor(matrix_t *A, int row, int col) {
  double det = 0;
  matrix_t minor;
  s21_create_matrix((*A).columns - 1, (*A).columns - 1, &minor);
  for (int i = 0; i < (*A).columns; i++) {
    int row_minor = i;
    if (row_minor > row) row_minor--;
    for (int j = 0; j < (*A).columns; j++) {
      int col_minor = j;
      if (col_minor > col) col_minor--;
      if (i != row && j != col)
        minor.matrix[row_minor][col_minor] = (*A).matrix[i][j];
    }
  }
  s21_determinant(&minor, &det);
  s21_remove_matrix(&minor);
  return det;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error = INCORRECT_MATRIX;
  if (s21_check_matrix(A)) {
    error = CALC_ERROR;
    if ((*A).columns == (*A).rows) {
      error = s21_create_matrix((*A).columns, (*A).rows, result);
      if (error == OK) {
        if ((*A).columns == 1) {
          (*result).matrix[0][0] = (*A).matrix[0][0];
        } else {
          for (int i = 0; i < (*A).columns; i++) {
            for (int j = 0; j < (*A).columns; j++) {
              (*result).matrix[i][j] =
                  s21_det_for_minor(A, i, j) * pow(-1, i + j);
            }
          }
        }
      }
    }
  }
  return error;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = INCORRECT_MATRIX;
  double det = 0;
  if (s21_check_matrix(A)) {
    error = CALC_ERROR;
    if (s21_determinant(A, &det) == 0) {
      if (det != 0) {
        if ((*A).columns == 1) {
          error = s21_create_matrix(1, 1, result);
          (*result).matrix[0][0] = 1.0 / (*A).matrix[0][0];
        } else {
          matrix_t matr_of_compl, inverse;
          s21_calc_complements(A, &matr_of_compl);
          s21_transpose(&matr_of_compl, &inverse);
          error = s21_mult_number(&inverse, 1.0 / det, result);
          s21_remove_matrix(&matr_of_compl);
          s21_remove_matrix(&inverse);
        }
      }
    }
  }
  return error;
}
