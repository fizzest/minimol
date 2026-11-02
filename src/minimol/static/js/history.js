// Undo/redo as a list of immutable models.

export class History {
  constructor(initial, limit = 200) {
    this.stack = [initial];
    this.pos = 0;
    this.limit = limit;
  }

  get current() {
    return this.stack[this.pos];
  }

  get canUndo() {
    return this.pos > 0;
  }

  get canRedo() {
    return this.pos < this.stack.length - 1;
  }

  /** Make `model` the current state. Returns false if it is already current. */
  push(model) {
    if (model === this.current) return false;
    this.stack.length = this.pos + 1; // drop the redo branch
    this.stack.push(model);
    if (this.stack.length > this.limit) this.stack.splice(0, this.stack.length - this.limit);
    this.pos = this.stack.length - 1;
    return true;
  }

  undo() {
    if (!this.canUndo) return null;
    this.pos -= 1;
    return this.current;
  }

  redo() {
    if (!this.canRedo) return null;
    this.pos += 1;
    return this.current;
  }
}
