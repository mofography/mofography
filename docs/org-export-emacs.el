;;; org-export-emacs.el — default dark HTML export for Org
;;; Add to ~/.emacs.d/config.el:
;;;   (load-file "~/repos/mofography/docs/org-export-emacs.el")

(defvar arun/org-html-css-file
  (expand-file-name "org-export-dark.css"
                    (expand-file-name "docs/" "~/repos/mofography/")))

(defun arun/org-html-styles ()
  (with-temp-buffer
    (insert-file-contents arun/org-html-css-file)
    (buffer-string)))

(defconst arun/org-html-copy-js
  "<script>
document.addEventListener('DOMContentLoaded',function(){
  document.querySelectorAll('.org-src-container').forEach(function(box){
    var pre=box.querySelector('pre');
    if(!pre||box.querySelector('.code-copy-btn'))return;
    box.classList.add('has-copy-btn');
    var btn=document.createElement('button');
    btn.type='button';btn.className='code-copy-btn';btn.textContent='Copy';
    btn.addEventListener('click',function(){
      var t=pre.innerText||pre.textContent||'';
      (navigator.clipboard&&navigator.clipboard.writeText
        ?navigator.clipboard.writeText(t)
        :Promise.reject(new Error('clipboard')))
      .then(function(){btn.textContent='Copied!';setTimeout(function(){btn.textContent='Copy'},1400);})
      .catch(function(){btn.textContent='Failed';setTimeout(function(){btn.textContent='Copy'},1400);});
    });
    box.insertBefore(btn,pre);
  });
});
</script>\n")

(setq org-html-head-include-default-style nil
      org-html-head-include-scripts nil
      org-html-validation-link nil
      org-babel-execute-on-export nil)

(setq org-html-head
      (concat "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
              "<style type=\"text/css\">\n"
              (arun/org-html-styles)
              "\n</style>\n"
              arun/org-html-copy-js))

(defun arun/org-html-strip-bad-headline-filters ()
  "Remove ALL custom headline/final-output filters (keep only org- builtins)."
  (dolist (sym '(org-export-filter-headline-functions
                 org-export-filter-final-output-functions))
    (when (boundp sym)
      (set sym
           (cl-remove-if
            (lambda (fn)
              (not (string-match-p "\\`org-" (symbol-name fn))))
            (symbol-value sym))))))

(defun arun/org-html-before-export (_backend)
  (arun/org-html-strip-bad-headline-filters))

(defun arun/org-html-sanitize-final-output (contents backend _info)
  "Repair HTML corrupted by § TOC permalink hooks."
  (unless (eq backend 'html)
    (cl-return-from arun/org-html-sanitize-final-output contents))
  ;; Strip § TOC permalink anchors
  (setq contents
        (replace-regexp-in-string
         "<a class=\"heading-permalink\"[^>]*>[^<]*</a>" "" contents))
  ;; Strip bare § TOC text that leaked into headings
  (setq contents
        (replace-regexp-in-string "§ TOC" "" contents))
  ;; Fix tag corruption
  (setq contents (replace-regexp-in-string "</hi>" "</h2>" contents))
  (setq contents (replace-regexp-in-string "<hiv " "<div " contents))
  (setq contents
        (replace-regexp-in-string "outline-co>tainer" "outline-container" contents))
  ;; Fix corrupted heading IDs: id="org<a class=..." → id="org..."
  (setq contents
        (replace-regexp-in-string
         "\\(id=\"\\)org[0-9a-f]*<a class=\"[^\"]*\"[^>]*>[^<]*</a>\""
         "\\1orgFIXED\""
         contents))
  (setq contents
        (replace-regexp-in-string
         "\\(id=\"\\)<a class=\"[^\"]*\"[^>]*>[^<]*</a>\""
         "\\1\""
         contents))
  ;; Drop empty / corrupted h2 left after permalink removal
  (setq contents
        (replace-regexp-in-string
         "<div id=\"outline-container-org\\([0-9a-f]+\\)\" class=\"outline-\\([0-9]+\\)\">\\s-*<div[^>]*>\\s-*<h2[^>]*>\\s-*</h2>"
         "<div id=\"outline-container-org\\1\" class=\"outline-\\2\">"
         contents))
  contents)
(arun/org-html-strip-bad-headline-filters)

(add-hook 'org-export-before-processing-hook #'arun/org-html-before-export)

(when (boundp 'org-export-filter-final-output-functions)
  (add-to-list 'org-export-filter-final-output-functions
               #'arun/org-html-sanitize-final-output))

(provide 'arun-org-export)
