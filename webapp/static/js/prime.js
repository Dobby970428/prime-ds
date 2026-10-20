/* ==========================================================================
   PRIME-DS — 프론트엔드 스크립트
   Streamlit 컴포넌트 래퍼 없이 Ketcher iframe과 직접 통신합니다.
   ========================================================================== */
(function () {
  "use strict";

  var CFG = window.PRIME || {};
  var KETCHER_URL = CFG.ketcherUrl || "https://ketcher.mireklzicar.com/";
  var STRUCT_URL = CFG.structureUrl || "/api/structure.png";
  var S = CFG.strings || {};

  /* ---------------- 유틸 ---------------- */

  function $(sel, root) { return (root || document).querySelector(sel); }
  function $$(sel, root) {
    return Array.prototype.slice.call((root || document).querySelectorAll(sel));
  }

  function toast(msg, kind, ms) {
    var wrap = $("#toasts");
    if (!wrap) return;
    var el = document.createElement("div");
    el.className = "toast" + (kind ? " toast--" + kind : "");
    el.textContent = msg;
    wrap.appendChild(el);
    setTimeout(function () {
      el.style.opacity = "0";
      el.style.transition = "opacity .25s ease";
      setTimeout(function () { el.remove(); }, 260);
    }, ms || 3800);
  }

  function overlay(on, title, text) {
    var el = $("#overlay");
    if (!el) return;
    if (title) $("#overlay-title").textContent = title;
    if (text) $("#overlay-text").textContent = text;
    el.classList.toggle("is-on", !!on);
  }

  window.PrimeUI = { toast: toast, overlay: overlay };

  /* ---------------- 사이드바 토글 ---------------- */

  var navtoggle = $("#navtoggle");
  if (navtoggle) {
    navtoggle.addEventListener("click", function () {
      var sb = $("#sidebar");
      if (sb) sb.classList.toggle("is-open");
    });
  }

  /* ==========================================================================
     Ketcher 구조 편집기
     --------------------------------------------------------------------------
     Ketcher iframe과의 통신 규약 (Streamlit 판에서 그대로 가져옴):
       보내기 : {type:"getSmiles"}  {type:"setMolecule", payload:"<smiles|molfile>"}
       받기   : {type:"init"}  {type:"smiles-update", smiles:"..."}
                {type:"smiles", payload:"..."}
     ========================================================================== */

  function Editor(root) {
    this.root = root;
    this.name = root.getAttribute("data-editor");
    this.frame = $('[data-role="frame"]', root);
    this.input = $('[data-role="input"]', root);
    this.status = $('[data-role="status"]', root);
    this.smilesEl = $('[data-role="smiles"]', root);
    this.preview = $('[data-role="preview"]', root);
    this.clearBtn = $('[data-role="clear"]', root);

    this.ready = false;
    this.pending = "";
    this.last = (this.input && this.input.value) || "";
    this.pollTimer = null;
    this.debounce = null;

    var self = this;

    this.frame.src = KETCHER_URL;

    if (this.clearBtn) {
      this.clearBtn.addEventListener("click", function () {
        self.post({ type: "setMolecule", payload: "" });
        self.apply("");
      });
    }

    /* ChemDraw 붙여넣기 — 캔버스에 포커스가 있을 때 Ctrl+V */
    root.addEventListener("paste", function (e) {
      var text = "";
      try {
        text = (e.clipboardData && e.clipboardData.getData("text/plain")) || "";
      } catch (err) { return; }
      if (!text || !text.trim()) return;
      var payload = text.trim();
      var isMol = payload.indexOf("M  END") >= 0 ||
                  payload.indexOf("V2000") >= 0 ||
                  payload.indexOf("V3000") >= 0;
      var isInchi = payload.indexOf("InChI=") === 0;
      var looksSmiles = /^\S{2,}$/.test(payload.split("\n")[0]);
      if (isMol || isInchi || looksSmiles) {
        e.preventDefault();
        self.load(payload);
        self.setStatus(S.pasting || "Applying paste…", "");
      }
    }, true);

    if (this.last) this.pending = this.last;
  }

  Editor.prototype.post = function (msg) {
    if (this.frame && this.frame.contentWindow) {
      this.frame.contentWindow.postMessage(msg, "*");
    }
  };

  Editor.prototype.load = function (payload) {
    if (this.ready) {
      this.post({ type: "setMolecule", payload: payload });
      var self = this;
      setTimeout(function () { self.post({ type: "getSmiles" }); }, 450);
    } else {
      this.pending = payload;
    }
  };

  Editor.prototype.setStatus = function (text, cls) {
    if (!this.status) return;
    this.status.textContent = text;
    this.status.className = "editor__status" + (cls ? " " + cls : "");
  };

  Editor.prototype.apply = function (smiles) {
    smiles = (smiles || "").trim();
    if (smiles === this.last) return;
    this.last = smiles;
    if (this.input) this.input.value = smiles;

    if (this.smilesEl) {
      this.smilesEl.textContent = smiles || S.placeholder || "Draw a structure";
    }
    if (smiles) {
      this.setStatus("✓ " + (S.saved || "Saved"), "is-ready");
    } else {
      this.setStatus(S.empty || "Empty", "is-empty");
    }
    if (this.preview) {
      if (smiles) {
        var img = new Image();
        img.src = STRUCT_URL + "?smiles=" + encodeURIComponent(smiles) + "&size=260x150";
        img.alt = S.noPreview || "preview";
        this.preview.innerHTML = "";
        this.preview.appendChild(img);
      } else {
        this.preview.innerHTML = '<span class="muted small">' + (S.noPreview || "No preview") + '</span>';
      }
    }
    /* 서버에 캔버스 상태를 조용히 저장 — 페이지를 옮겨도 유지되도록 */
    this.sync(smiles);
    this.root.dispatchEvent(new CustomEvent("prime:smiles", {
      bubbles: true, detail: { name: this.name, smiles: smiles }
    }));
  };

  Editor.prototype.sync = function (smiles) {
    clearTimeout(this.debounce);
    var self = this;
    this.debounce = setTimeout(function () {
      fetch("/api/canvas", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ name: self.name, smiles: smiles })
      }).catch(function () { /* silent */ });
    }, 500);
  };

  Editor.prototype.handle = function (data) {
    if (!data || typeof data !== "object") return;

    if (data.type === "init" || data.eventType === "init") {
      this.ready = true;
      var self = this;
      if (!this.pollTimer) {
        this.pollTimer = setInterval(function () {
          if (self.ready) self.post({ type: "getSmiles" });
        }, 2200);
      }
      if (this.pending) {
        this.post({ type: "setMolecule", payload: this.pending });
        this.pending = "";
      }
      this.post({ type: "getSmiles" });
      return;
    }

    var smi = null;
    if (data.type === "smiles-update" && typeof data.smiles === "string") {
      smi = data.smiles;
    } else if (data.type === "smiles" && typeof data.payload === "string") {
      smi = data.payload;
    }
    if (smi !== null) this.apply(smi);
  };

  var editors = [];
  $$("[data-editor]").forEach(function (el) { editors.push(new Editor(el)); });

  if (editors.length) {
    window.addEventListener("message", function (event) {
      /* 어느 iframe에서 온 메시지인지 source로 구분 */
      for (var i = 0; i < editors.length; i++) {
        if (editors[i].frame && editors[i].frame.contentWindow === event.source) {
          editors[i].handle(event.data);
          return;
        }
      }
    });
  }

  window.PrimeEditors = {
    get: function (name) {
      for (var i = 0; i < editors.length; i++) {
        if (editors[i].name === name) return editors[i];
      }
      return null;
    },
    smiles: function (name) {
      var ed = this.get(name);
      return ed ? ed.last : "";
    },
    load: function (name, payload) {
      var ed = this.get(name);
      if (ed) ed.load(payload);
    }
  };

  /* ==========================================================================
     폼 — 제출 시 오버레이
     ========================================================================== */

  $$("form[data-busy]").forEach(function (form) {
    form.addEventListener("submit", function () {
      var title = form.getAttribute("data-busy") || "Working…";
      var text = form.getAttribute("data-busy-text") || "";
      overlay(true, title, text);
    });
  });

  /* 제출 전 필수 캔버스 확인 */
  $$("form[data-require-canvas]").forEach(function (form) {
    form.addEventListener("submit", function (e) {
      var names = form.getAttribute("data-require-canvas").split(",");
      var missing = [];
      names.forEach(function (n) {
        n = n.trim();
        if (!n) return;
        var input = form.querySelector('input[name="' + n + '"]');
        if (!input || !input.value.trim()) missing.push(n);
      });
      if (missing.length) {
        e.preventDefault();
        overlay(false);
        toast((S.drawFirst || "Draw a structure first —") + " " + missing.join(", "), "warn");
      }
    });
  });

  /* ==========================================================================
     range 슬라이더 값 표시
     ========================================================================== */

  $$('input[type="range"][data-output]').forEach(function (input) {
    var out = document.getElementById(input.getAttribute("data-output"));
    if (!out) return;
    var suffix = input.getAttribute("data-suffix") || "";
    function render() { out.textContent = input.value + suffix; }
    input.addEventListener("input", render);
    render();
  });

  /* ==========================================================================
     의존 필드 — 특정 select 값에 따라 다른 입력을 비활성화
     ========================================================================== */

  $$("[data-enable-when]").forEach(function (target) {
    var spec = target.getAttribute("data-enable-when").split("!=");
    var srcName = spec[0].trim();
    var forbidden = (spec[1] || "").trim();
    var src = document.querySelector('[name="' + srcName + '"]');
    if (!src) return;
    function render() {
      var off = src.value === forbidden;
      target.disabled = off;
      target.closest(".field") &&
        target.closest(".field").classList.toggle("is-off", off);
    }
    src.addEventListener("change", render);
    render();
  });

  /* ==========================================================================
     탭
     ========================================================================== */

  $$("[data-tabs]").forEach(function (group) {
    var tabs = $$(".tab", group);
    tabs.forEach(function (tab) {
      tab.addEventListener("click", function () {
        var target = tab.getAttribute("data-target");
        tabs.forEach(function (t) { t.classList.remove("is-active"); });
        tab.classList.add("is-active");
        var scope = group.parentElement || document;
        $$(".tabpanel", scope).forEach(function (p) {
          p.classList.toggle("is-active", p.id === target);
        });
      });
    });
  });

  /* ==========================================================================
     부형제 선택 — 펼치기 · 접기 · 전체 선택
     ========================================================================== */

  $$("[data-exc-toggle]").forEach(function (btn) {
    btn.addEventListener("click", function () {
      var open = btn.getAttribute("data-exc-toggle") === "open";
      $$(".exc-group").forEach(function (d) { d.open = open; });
    });
  });

  $$("[data-exc-selectall]").forEach(function (btn) {
    btn.addEventListener("click", function () {
      var boxes = $$('#exc-pick input[type="checkbox"][name="key"]:not(:disabled)');
      var allOn = boxes.every(function (b) { return b.checked; });
      boxes.forEach(function (b) { b.checked = !allOn; });
    });
  });

  /* ==========================================================================
     클립보드 복사
     ========================================================================== */

  document.addEventListener("click", function (e) {
    var btn = e.target.closest("[data-copy]");
    if (!btn) return;
    var text = btn.getAttribute("data-copy");
    if (!text) {
      var sel = btn.getAttribute("data-copy-from");
      var el = sel && document.querySelector(sel);
      text = el ? (el.value || el.textContent) : "";
    }
    if (!text) return;
    navigator.clipboard.writeText(text).then(
      function () { toast(S.copied || "Copied", "ok", 1600); },
      function () { toast(S.copyFailed || "Copy failed", "danger"); }
    );
  });

  /* ==========================================================================
     페이지 이탈 시 오버레이 정리 (뒤로가기 대응)
     ========================================================================== */

  window.addEventListener("pageshow", function () { overlay(false); });
})();
